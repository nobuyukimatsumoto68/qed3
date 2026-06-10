#!/usr/bin/env bash
# tmp_AB_claude.sh -- free L1 methods A and B for the CFT comparison.
#   A: exact K (basis loop) + lattice P.  RECOMPILE jj_contract_deter (antiperiodic conn_shift fix) and
#      re-contract on the CACHED K (data_free_Kdense_L1/K.0.h5) + cached propagator -> corr_deter_exact_L1.
#      (No K rebuild: K(0) is correct; the BC fix is only in the time-shift.)
#   B: local op (BARE-GAMMA redesign) + lattice P.  RECOMPILE jj_local_deter (the on-disk .o is the old
#      Omega/build_W version) and run -> corr_deter_local_L1 (overwrites the stale file).
# Both GPU.  Pick the GPU below.  Run:  bash tmp_AB_claude.sh 2>&1 | tee AB_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES=0          # <-- set to the free GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

D=data_free_vmRe0.000000vmIm0.000000
[ -f data_free_Kdense_L1/K.0.h5 ]   || { echo "missing cached exact K (data_free_Kdense_L1/K.0.h5)"; exit 1; }
[ -f $D/prop_deter_L1/Dinv.0.h5 ]   || { echo "missing cached propagator $D/prop_deter_L1/Dinv.0.h5"; exit 1; }

echo "===== A: compile jj_contract_deter_L1.o (antiperiodic conn_shift) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=1 jj_contract_deter_claude.cu -o jj_contract_deter_L1.o \
  || { echo "BUILD FAILED (contract)"; exit 1; }
echo "===== A: run exact contraction (free, L1, massless) on CACHED K ====="
rm -f $D/corr_deter_exact_L1/corr.0.h5*
./jj_contract_deter_L1.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
  || { echo "CONTRACT FAILED"; exit 1; }

echo "===== B: compile jj_local_deter_L1.o (bare-gamma redesign) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=1 jj_local_deter_claude.cu -o jj_local_deter_L1.o \
  || { echo "BUILD FAILED (local)"; exit 1; }
echo "===== B: run local op (free, L1, massless) ====="
rm -f $D/corr_deter_local_L1/corr.0.h5*
./jj_local_deter_L1.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
  || { echo "LOCAL FAILED"; exit 1; }

echo "===== quick 3-way tp shape + 2-way sp shape (Python) ====="
OMP_NUM_THREADS=1 python3 - << 'PY'
import h5py, numpy as np
D='data_free_vmRe0.000000vmIm0.000000'
def vpp(path, proj, key='h0/t0_0/'+'%s/Vpp'):
    f=h5py.File(path,'r'); k=('h0/t0_0/%s/Vpp')%proj
    return f[k+'/real'][:]+1j*f[k+'/imag'][:]
def gt_avg(path):
    f=h5py.File(path,'r'); return f['h0/tp/Gt_avg/real'][:]+1j*f['h0/tp/Gt_avg/imag'][:]
ds=[1,2,3,4,6,8,12,16,24,32]
ex=D+'/corr_deter_exact_L1/corr.0.h5'
lo=D+'/corr_deter_local_L1/corr.0.h5'
co=D+'/corr_deter_commut_L1/corr.0.h5'
import os
series={}
if os.path.exists(ex): series['exact_tp']=vpp(ex,'tp').real; series['exact_sp']=vpp(ex,'sp').real
if os.path.exists(lo): series['local_tp']=vpp(lo,'tp').real; series['local_sp']=vpp(lo,'sp').real
if os.path.exists(co): series['commut_tp']=gt_avg(co).real
def shape(g):   # subtract large-dt plateau, normalize to dt=1
    c=g-np.mean(g[40:80]); return c/c[1]
print(" dt      :", " ".join("%7d"%d for d in ds))
for name,g in series.items():
    s=shape(np.asarray(g))
    print(" %-9s:"%name, " ".join("%+7.3f"%s[d] for d in ds))
print("\n(tp shapes exact/local/commut should overlay; sp shapes exact/local should overlay)")
PY
echo "===== DONE ====="
