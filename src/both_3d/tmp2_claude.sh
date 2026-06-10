#!/usr/bin/env bash
# tmp2_claude.sh -- LOCAL-current validation on GPU 1 (GPU 0 is running the exact test).
# Builds ONLY jj_local_deter (does NOT recompile the propagator binary GPU 0 is using), reuses the
# existing free propagator data_free_.../prop_deter_L1/Dinv.0.h5, runs the local current (free, L=1,
# massless), and compares LOCAL vs EXACT (if corr_deter_exact_L1 is ready).
# Run:  bash tmp2_claude.sh 2>&1 | tee local_exact_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

PROP=data_free_vmRe0.000000vmIm0.000000/prop_deter_L1/Dinv.0.h5
[ -f "$PROP" ] || { echo "missing $PROP (the C1 propagator) -- run that first"; exit 1; }

echo "===== compile jj_local_deter_L1.o (only) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=1 jj_local_deter_claude.cu -o jj_local_deter_L1.o \
  || { echo "BUILD FAILED"; exit 1; }

echo "===== run LOCAL current (free, L=1, massless) on GPU 1 ====="
# validation: always recompute (remove the prior output so a code change is actually picked up;
# the 'complete' sentinel would otherwise make the run skip).
rm -f data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L1/corr.0.h5*
CUDA_VISIBLE_DEVICES=1 ./jj_local_deter_L1.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
  || { echo "RUN FAILED"; exit 1; }

echo "===== compare LOCAL vs EXACT  (G_s/G_t ratio; CFT: -(D-1)=-2) ====="
LOC=data_free_vmRe0.000000vmIm0.000000/corr_deter_local_L1/corr.0.h5
EXA=data_free_vmRe0.000000vmIm0.000000/corr_deter_exact_L1/corr.0.h5
python3 - "$LOC" "$EXA" << 'PY'
import h5py, sys, numpy as np, os
def load(p):
    f=h5py.File(p,'r')
    g=lambda L: f[L+'/real'][:]+1j*f[L+'/imag'][:]
    return g('h0/t0_0/tp/Vpp').real, g('h0/t0_0/sp/Vpp').real
ds=[1,2,3,4,6,8,12,16]
lt,ls=load(sys.argv[1]); print("LOCAL  loaded")
print(" dt   :"," ".join("%7d"%d for d in ds))
print(" Gt(loc):"," ".join("%+7.1e"%lt[d] for d in ds))
print(" Gs(loc):"," ".join("%+7.1e"%ls[d] for d in ds))
with np.errstate(divide='ignore',invalid='ignore'): r=ls/lt
print(" Gs/Gt  :"," ".join("%+7.3f"%r[d] for d in ds),"  (CFT -2)")
if os.path.exists(sys.argv[2]):
    et,es=load(sys.argv[2]); print("\nEXACT loaded -- local/exact ratio (should be ~const if local is a good proxy):")
    with np.errstate(divide='ignore',invalid='ignore'):
        print(" Gt loc/exa:"," ".join("%+7.3f"%(lt[d]/et[d]) for d in ds))
        print(" Gs loc/exa:"," ".join("%+7.3f"%(ls[d]/es[d]) for d in ds))
else:
    print("\n(exact corr_deter_exact_L1 not ready yet -- rerun this compare once GPU0's run finishes)")
PY
echo "===== DONE ====="
