#!/usr/bin/env bash
# Build + smoke-test the exact propagator builder (C1), L=1, free-field massless.
# Run:  bash tmp_claude.sh 2>&1 | tee prop_exact_build_claude.log
set -u
LOG_TAG(){ echo; echo "===== $* ====="; }

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null

LOG_TAG "BUILD jj_propagator_exact_claude.o (L=1, default N_REFINE_CLI=1)"
make jj_propagator_exact_claude.o || { echo "BUILD FAILED"; exit 1; }

LOG_TAG "RUN free-field massless (U=1) on GPU0"
# free field = no --ens-dir; massless = mass 0; writes data_free_vmRe0.000000vmIm0.000000/prop_exact_L1/Dinv.0.h5
CUDA_VISIBLE_DEVICES=0 ./jj_propagator_exact_claude.o --mass-re 0.0 --mass-im 0.0 || { echo "RUN FAILED"; exit 1; }

LOG_TAG "OUTPUT check"
OUT=data_free_vmRe0.000000vmIm0.000000/prop_exact_L1/Dinv.0.h5
ls -l "$OUT" 2>/dev/null && python3 - "$OUT" << 'PY'
import h5py, sys, numpy as np
f=h5py.File(sys.argv[1],'r')
print("keys:", list(f.keys()))
N=f['N'][0]; print("N=",N,"parity=",f['parity'][0],"complete=",'complete' in f)
re=np.array(f['Dm_inv/real']); im=np.array(f['Dm_inv/imag'])
print("Dm_inv len=",re.size," expect N^2=",N*N," match=",re.size==N*N)
P=(re+1j*im).reshape(N,N)
print("Dm_inv[0,0]=",P[0,0]," |P| frob=",np.linalg.norm(P))
PY
LOG_TAG "DONE"
