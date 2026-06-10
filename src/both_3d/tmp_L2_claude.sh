#!/usr/bin/env bash
# tmp_L2_claude.sh -- free L=2 (N=10752): build the dense propagator (Dov + P) then method B (local op).
# Reduces cutoff effects vs L=1.  GPU 0 (Titan V, 12 GB -> L2 dense matrices 1.85 GB each fit; L4 does NOT).
# Run:  bash tmp_L2_claude.sh 2>&1 | tee L2_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES=0

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "===== compile L=2 propagator + local (-DN_REFINE_CLI=2) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=2 jj_propagator_deter_claude.cu -o jj_propagator_deter_L2.o \
  || { echo "BUILD FAILED (propagator)"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=2 jj_local_deter_claude.cu -o jj_local_deter_L2.o \
  || { echo "BUILD FAILED (local)"; exit 1; }

echo "===== [1] build dense propagator (free, L2, massless) -- Dov + P + LU ====="
./jj_propagator_deter_L2.o --mass-re 0.0 --mass-im 0.0 \
  || { echo "PROPAGATOR FAILED"; exit 1; }

echo "===== [2] method B: local op + lattice P (free, L2) ====="
rm -f data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L2/corr.0.h5*
./jj_local_deter_L2.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
  || { echo "LOCAL FAILED"; exit 1; }

echo "===== B(L2) ratio Gs/Gt (connected) ====="
OMP_NUM_THREADS=1 python3 - << 'PY'
import h5py, numpy as np
p='data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L2/corr.0.h5'
f=h5py.File(p,'r')
def v(proj): k='h0/t0_0/%s/Vpp'%proj; return f[k+'/real'][:]+1j*f[k+'/imag'][:]
def conn(g): return (g-np.mean(g[40:80])).real
gt,gs=conn(v('tp')),conn(v('sp'))
ds=[1,2,3,4,6,8]
with np.errstate(divide='ignore',invalid='ignore'): r=gs/gt
print(" dt   :", " ".join("%8d"%d for d in ds))
print(" Gs/Gt:", " ".join("%+8.3f"%r[d] for d in ds), "  (CFT -2)")
PY
echo "===== DONE ====="
