#!/bin/sh
#SBATCH --account=qed3.lq2_gpu
#SBATCH --qos=normal
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --mem=192G    # host (CPU) RAM, NOT GPU mem (GPU is set by --gpus). Covers stage 2's host-dense
#                     # contraction (~110 GB peak; README 128-256 GB). Explicit -> deterministic footprint.
#SBATCH --time=12:00:00   # stage 1 is ~10-30 min, but stage 2/3 host-dense tr(A(t0)A(t)) at L=4 (O(N^2 Nt),
#                         # + reading the 55 GB P) can take hours; qos=normal allows a 1-day wall, so be generous.
#SBATCH --job-name=L4_a100_freeprop
#SBATCH --output=slurm_%x_%j.out

# L4_a100_batch_claude.sh -- SLURM batch: free L=4 current-correlator on an A100-80.
# See L4_a100_README_claude.md for goal, memory budget, and what to return.
#   Stage 1 (A100): dense overlap propagator  L4_a100_prop_4.o     -> data_free_.../prop_deter_L4/Dinv.0.h5  (~55 GB)
#   Stage 2 (host): B = local op + lattice P   jj_local_deter_4.o   -> data_free_.../corr_deter_local_L4/corr.0.h5
#   Stage 3 (host): D = local op + continuum G jj_local_deter_4.o   -> data_free_.../corr_deter_local_cont_L4/corr.0.h5
#                   (same binary, --prop-file cont_prop_L4/Dinv.0.h5 --out-tag cont; auto-skips if that file is absent)
# Binaries are built by L4_a100_make_submit_claude.sh (which submits this script).
# A100-80 needed: stage 1 holds D_ov (27.5 GB) + its inverse (27.5 GB) + LU workspace \approx 56 GB on the GPU.
# Host RAM: stage 2 host-dense contraction needs ~110 GB (P + A(t0) + A(t)); --mem=192G below covers it.

set -u

# print worker node + environment
hostname
source /lustre2/qed3/env.sh
date
nvidia-smi

cd /project/qed3/qed3/src/both_3d || { echo "cd to src/both_3d FAILED"; exit 1; }

L=4

echo "===== [1] dense propagator (free, L=${L}, massless) -- Dov + P + LU (A100, ~56 GB GPU) ====="
# expect: "[check] || D_m . D_m^{-1} - I ||_F = ~1e-13"  ->  data_free_.../prop_deter_L${L}/Dinv.0.h5 (~55 GB)
time ./L4_a100_prop_${L}.o --mass-re 0.0 --mass-im 0.0 || { echo "PROPAGATOR FAILED"; exit 1; }

echo "===== [2] method B: local op + lattice P (free, L=${L}) -- host-dense (~110 GB RAM) ====="
# NO rm here: jj_local_deter self-skips a COMPLETE corr.0.h5 (its own 'complete' sentinel, src:301-303) and
# truncates any stale .tmp on open (src:309).  A pre-emptive 'rm corr.0.h5*' would DESTROY a finished result
# on resubmit (it did once -- wiped a 109-min run).  Let the binary's guard do its job.
time ./jj_local_deter_${L}.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 || { echo "LOCAL FAILED (host RAM? -> GPU port, see README)"; exit 1; }

echo "===== B(L${L}) ratio Gs/Gt (connected) -- CFT target -2 ====="
OMP_NUM_THREADS=1 python3 - << 'PY'
import h5py, numpy as np
p='data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L4/corr.0.h5'
f=h5py.File(p,'r')
def v(proj): k='h0/t0_0/%s/Vpp'%proj; return f[k+'/real'][:]+1j*f[k+'/imag'][:]
def conn(g): return (g-np.mean(g[40:80])).real
gt,gs=conn(v('tp')),conn(v('sp'))
with np.errstate(divide='ignore',invalid='ignore'): r=gs/gt
print(" dt   :", " ".join("%8d"%d for d in [1,2,3,4,6,8]))
print(" Gs/Gt:", " ".join("%+8.3f"%r[d] for d in [1,2,3,4,6,8]), "  (CFT -2)")
PY

echo "===== [3] method D: local op + CONTINUUM G (free, L=${L}) -- THE CFT verification ====="
# Same jj_local_deter binary, but read P from the continuum propagator cont_prop_L4/Dinv.0.h5 (same Dinv
# schema) via --prop-file, writing corr_deter_local_cont_L${L}/.  Auto-skips if that file is absent (README).
CONT_PROP=cont_prop_L${L}/Dinv.0.h5
if [ -f "${CONT_PROP}" ]; then
  # NO rm (see stage 2): the binary self-skips a complete corr.0.h5 and truncates stale .tmp itself.
  time ./jj_local_deter_${L}.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
       --prop-file "${CONT_PROP}" --out-tag cont || { echo "D (continuum) FAILED"; exit 1; }
  echo "===== D(L${L}) ratio Gs/Gt (connected) -- CFT target -2 ====="
  OMP_NUM_THREADS=1 python3 - << 'PY'
import h5py, numpy as np
p='data_free_vmRe0.000000vmIm0.000000/corr_deter_local_cont_L4/corr.0.h5'
f=h5py.File(p,'r')
def v(proj): k='h0/t0_0/%s/Vpp'%proj; return f[k+'/real'][:]+1j*f[k+'/imag'][:]
def conn(g): return (g-np.mean(g[40:80])).real
gt,gs=conn(v('tp')),conn(v('sp'))
with np.errstate(divide='ignore',invalid='ignore'): r=gs/gt
print(" dt   :", " ".join("%8d"%d for d in [1,2,3,4,6,8]))
print(" Gs/Gt:", " ".join("%+8.3f"%r[d] for d in [1,2,3,4,6,8]), "  (CFT -2)")
PY
else
  echo "# SKIP stage 3: ${CONT_PROP} not present (transfer the 27.5 GB continuum propagator to run method D)."
fi

echo "===== DONE -- return data_free_vmRe0.000000vmIm0.000000/{corr_deter_local_L${L},corr_deter_local_cont_L${L}}/corr.0.h5 ====="
exit
