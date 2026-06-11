#!/usr/bin/env bash
# tmp_diag_claude.sh -- rebuild jj_local_deter (now INSERTION-DIAGONAL, Eq. 4.29 sp / 4.32 tp) and
# rerun free L=1 massless -> corr_deter_local_L1 (overwrites).  Then compare to the exact stochastic K:
# local sp should now be POSITIVE (Eq. 4.28, matches exact), local tp NEGATIVE in the IR (Eq. 4.31).
# Run:  bash tmp_diag_claude.sh 2>&1 | tee diag_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES=0          # <-- set to a free GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

D=data_free_vmRe0.000000vmIm0.000000
[ -f $D/prop_deter_L1/Dinv.0.h5 ] || { echo "missing cached propagator $D/prop_deter_L1/Dinv.0.h5"; exit 1; }

echo "===== compile jj_local_deter_L1.o (insertion-diagonal) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=1 jj_local_deter_claude.cu -o jj_local_deter_L1.o \
  || { echo "BUILD FAILED"; exit 1; }
echo "===== run local op DIAGONAL (free, L1, massless) ====="
rm -f $D/corr_deter_local_L1/corr.0.h5*
./jj_local_deter_L1.o --mass-re 0.0 --mass-im 0.0 --n-t0 2 \
  || { echo "LOCAL FAILED"; exit 1; }

echo "===== sign check: exact (diagonal stochastic) vs local (diagonal) ====="
OMP_NUM_THREADS=1 python3 - << 'PY'
import h5py, numpy as np, os
D='data_free_vmRe0.000000vmIm0.000000'
def vpp_groups(path, proj, origin='t0_0'):   # hit-averaged over h0/,h1/,...
    out=[]
    with h5py.File(path,'r') as f:
        hs=sorted([k for k in f if k.startswith('h') and k[1:].isdigit()], key=lambda s:int(s[1:]))
        for hg in hs:
            k=f'{hg}/{origin}/{proj}/Vpp'
            if (k+'/real') in f: out.append(f[k+'/real'][()]+1j*f[k+'/imag'][()])
    return np.array(out).mean(0)
exact=D+'/corr_nt02_nhits64/corr.0.h5'
local=D+'/corr_deter_local_L1/corr.0.h5'
print(" dt          :", " ".join("%9d"%d for d in range(1,7)))
for name,path in [('exact',exact),('local',local)]:
    if not os.path.exists(path): print(" %-6s: MISSING"%name); continue
    for proj in ['tp','sp']:
        g=vpp_groups(path,proj).real
        print(" %-6s %s:"%(name,proj), " ".join("%+9.2e"%g[d] for d in range(1,7)))
print("\nEXPECT: local sp now POSITIVE (was negative) and same sign as exact sp;")
print("        local tp negative in the IR (Eq. 4.31).  Magnitudes differ (bare W vs dressed K).")
PY
echo "===== DONE ====="
