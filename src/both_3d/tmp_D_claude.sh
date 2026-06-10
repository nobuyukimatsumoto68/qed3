#!/usr/bin/env bash
# tmp_D_claude.sh -- method D (CFT verification): local op (bare-gamma) + CONTINUUM propagator G at L1.
# Recompiles jj_local_deter (path wiring for --prop-file/--out-tag) and runs it reading the continuum
# free propagator cont_prop_L1/Dinv.0.h5 -> corr_deter_local_cont_L1.  Then checks the parameter-free CFT
# ratio Gs/Gt = -(D-1) = -2 and the continuum-vs-lattice-local shapes.
# Run:  bash tmp_D_claude.sh 2>&1 | tee D_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES=0          # <-- set to the free GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

[ -f cont_prop_L1/Dinv.0.h5 ] || { echo "missing continuum prop cont_prop_L1/Dinv.0.h5"; exit 1; }

echo "===== compile jj_local_deter_L1.o (with --prop-file/--out-tag wiring) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=1 jj_local_deter_claude.cu -o jj_local_deter_L1.o \
  || { echo "BUILD FAILED"; exit 1; }

echo "===== D: run local op + CONTINUUM G (free, L1, massless) ====="
rm -f data_free_vmRe0.000000vmIm0.000000/corr_deter_local_cont_L1/corr.0.h5*
./jj_local_deter_L1.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
  --prop-file cont_prop_L1/Dinv.0.h5 --out-tag cont \
  || { echo "RUN FAILED"; exit 1; }

echo "===== CFT ratio Gs/Gt = -(D-1) = -2  +  continuum-vs-lattice shapes ====="
OMP_NUM_THREADS=1 python3 - << 'PY'
import h5py, numpy as np, os
D='data_free_vmRe0.000000vmIm0.000000'
def vpp(path, proj):
    f=h5py.File(path,'r'); k='h0/t0_0/%s/Vpp'%proj
    return (f[k+'/real'][:]+1j*f[k+'/imag'][:]).real
cont=D+'/corr_deter_local_cont_L1/corr.0.h5'
latt=D+'/corr_deter_local_L1/corr.0.h5'
gt=vpp(cont,'tp'); gs=vpp(cont,'sp')
def conn(g): return g-np.mean(g[40:80])   # subtract large-dt plateau
ct,cs=conn(gt),conn(gs)
ds=[1,2,3,4,6,8,12,16]
print(" dt        :", " ".join("%8d"%d for d in ds))
print(" Gt_conn   :", " ".join("%+8.2e"%ct[d] for d in ds))
print(" Gs_conn   :", " ".join("%+8.2e"%cs[d] for d in ds))
with np.errstate(divide='ignore',invalid='ignore'):
    r=cs/ct
print(" Gs/Gt     :", " ".join("%+8.3f"%r[d] for d in ds), "   (CFT: -(D-1) = -2)")
if os.path.exists(latt):
    lt=conn(vpp(latt,'tp')); ls=conn(vpp(latt,'sp'))
    print("\n continuum vs lattice-local (shapes, normalized to dt=1):")
    print(" tp cont   :", " ".join("%+7.3f"%(ct[d]/ct[1]) for d in ds))
    print(" tp latt   :", " ".join("%+7.3f"%(lt[d]/lt[1]) for d in ds))
    print(" sp cont   :", " ".join("%+7.3f"%(cs[d]/cs[1]) for d in ds))
    print(" sp latt   :", " ".join("%+7.3f"%(ls[d]/ls[1]) for d in ds))
print("\n NOTE: L=1 is only 12 sites (coarse); the ratio -> -2 and CFT shape sharpen at L=2,4.")
PY
echo "===== DONE ====="
