#!/usr/bin/env bash
# L4_a100_run_claude.sh -- build + run the free L=4 DETER current-correlator computation on an A100-80.
# Stage 1 (A100): dense propagator   L4_a100_prop_claude.cu   -> data_free_.../prop_deter_L4/Dinv.0.h5
# Stage 2 (host): B = local + lattice P   jj_local_deter_claude.cu -> .../corr_deter_local_L4/corr.0.h5
# Stage 3 (host): D = local + continuum G  (same binary, --prop-file cont_prop_L4) -> .../corr_deter_local_cont_L4/
# See L4_a100_README_claude.md for prerequisites, memory, and what to return.
# Run:  bash L4_a100_run_claude.sh 2>&1 | tee L4_a100_claude.log
set -u
cd "$(dirname "$0")" || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES=0                       # <-- the A100

NVCC=nvcc
# NOTE: -arch=sm_80 for A100 (the repo default is sm_70 for Titan V).
NVCCFLAGS="-arch=sm_80 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
# ADJUST INCLUDES/LDFLAGS to this server's Eigen / HighFive / HDF5 locations.

L=4
echo "===== compile (N_REFINE_CLI=${L}, sm_80) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} L4_a100_prop_claude.cu   -o L4_a100_prop_${L}.o    || { echo "BUILD FAILED (prop)";  exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_local_deter_claude.cu -o jj_local_deter_${L}.o  || { echo "BUILD FAILED (local)"; exit 1; }

echo "===== [1] dense propagator (free, L4, massless) -- Dov + P + LU (A100, ~56 GB GPU) ====="
# expect: "[check] || D_m . D_m^{-1} - I ||_F = ~1e-13"  ->  data_free_.../prop_deter_L4/Dinv.0.h5 (~55 GB)
./L4_a100_prop_${L}.o --mass-re 0.0 --mass-im 0.0 || { echo "PROPAGATOR FAILED"; exit 1; }

echo "===== [2] method B: local op + lattice P (free, L4) -- host-dense (~110 GB RAM) ====="
# if host RAM is insufficient, this stage must be GPU-ported (see README); the propagator (stage 1) is the
# A100-critical part and is complete by here.
rm -f data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L4/corr.0.h5*
./jj_local_deter_${L}.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 || { echo "LOCAL FAILED (host RAM? -> GPU port)"; exit 1; }

echo "===== B(L4) ratio Gs/Gt (connected) -- CFT target -2 ====="
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
echo "===== [3] method D: local op + CONTINUUM G (free, L4) -- host-dense (~110 GB RAM) ====="
# needs cont_prop_L4/Dinv.0.h5 (27.5 GB) on this server. Same binary, propagator overridden.
if [ -f cont_prop_L4/Dinv.0.h5 ]; then
  rm -f data_free_vmRe0.000000vmIm0.000000/corr_deter_local_cont_L4/corr.0.h5*
  ./jj_local_deter_${L}.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
    --prop-file cont_prop_L4/Dinv.0.h5 --out-tag cont || { echo "D FAILED (host RAM? -> GPU port)"; exit 1; }
else
  echo "  SKIP stage 3: cont_prop_L4/Dinv.0.h5 not present on this server (transfer it, 27.5 GB)."
fi

echo "===== DONE -- return corr_deter_local_L4/corr.0.h5  (and corr_deter_local_cont_L4/corr.0.h5 if stage 3 ran) ====="
