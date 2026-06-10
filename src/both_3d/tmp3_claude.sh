#!/usr/bin/env bash
# tmp3_claude.sh -- DEBUG the broken C1 propagator on GPU 1: dump dense D_ov and test whether D_ov
# itself is time-translation invariant.  If D_ov is invariant but P=D_ov^{-1} is not => cuSOLVER LU
# bug; if D_ov is NOT invariant => the dense from_cpu build is wrong.
# Compiles to a DISTINCT binary (jj_propagator_deter_dbg.o) so it can't clobber the exact run's
# jj_propagator_deter_L1.o.  Regenerates the free Dinv.0.h5 (now with a Dov dump).
# Run:  bash tmp3_claude.sh 2>&1 | tee prop_dov_debug_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "===== compile jj_propagator_deter_dbg.o (with D_ov dump) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=1 jj_propagator_deter_claude.cu -o jj_propagator_deter_dbg.o \
  || { echo "BUILD FAILED"; exit 1; }

echo "===== regenerate free Dinv.0.h5 (rm stale; recompute with dump) on GPU 1 ====="
rm -f data_free_vmRe0.000000vmIm0.000000/prop_deter_L1/Dinv.0.h5*
CUDA_VISIBLE_DEVICES=1 ./jj_propagator_deter_dbg.o --mass-re 0.0 --mass-im 0.0 \
  || { echo "RUN FAILED"; exit 1; }

echo "===== translation-invariance test: D_ov vs D_ov^{-1} ====="
python3 - << 'PY'
import h5py, numpy as np
f=h5py.File('data_free_vmRe0.000000vmIm0.000000/prop_deter_L1/Dinv.0.h5','r')
N=int(f['N'][0]); Nt=128; Nx=N//Nt
def get(key): return (f[key+'/real'][:]+1j*f[key+'/imag'][:]).reshape(N,N)
def bn(M,t,dt): tp=(t+dt)%Nt; return np.linalg.norm(M[t*Nx:(t+1)*Nx,tp*Nx:(tp+1)*Nx])
for key in ['Dov','Dm_inv']:
    if key+'/real' not in f: print(f"  {key}: not present"); continue
    M=get(key); print(f"=== {key} ===")
    for dt in [0,1,2,4]:
        vals=[bn(M,t,dt) for t in range(8)]
        print(f"  dt={dt}: ||{key}[t,t+{dt}]|| t=0..7: "+" ".join("%.2e"%v for v in vals)
              +"   (constant over t => translation invariant)")
PY
echo "===== DONE ====="
