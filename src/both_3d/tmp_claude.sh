#!/usr/bin/env bash
# Drive the EXACT pipeline (free, L=1, massless) via the orchestrator, then check Eq.4.28.
# Stage 1b (K-build) is the slow part (~30-50 min, N x n_insertions multishift applies);
# the propagator (stage 1a) already exists from the earlier C1 run and will be skipped.
# Run:  bash tmp_claude.sh 2>&1 | tee contract_exact_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

echo "===== run_exact_claude.sh  (free, L=1, mass 0:0) ====="
bash run_exact_claude.sh --L 1 --masses "0:0" --free --gpu 0 --n-t0 2 || { echo "PIPELINE FAILED"; exit 1; }

echo
echo "===== Eq.4.28 ratio check  (CFT: G_s/G_t = -(D-1) = -2 at every dt) ====="
OUT=data_free_vmRe0.000000vmIm0.000000/corr_exact_L1/corr.0.h5
ls -l "$OUT" 2>/dev/null && python3 - "$OUT" << 'PY'
import h5py, sys, numpy as np
f=h5py.File(sys.argv[1],'r')
def g(leaf): return f[leaf+'/real'][:]+1j*f[leaf+'/imag'][:]
tp=g('h0/t0_0/tp/Vpp').real; sp=g('h0/t0_0/sp/Vpp').real
ds=[1,2,3,4,6,8,12,16]
print("complete=",'complete' in f," n_t0=",f['n_t0'][0])
print(" dt   :", " ".join("%7d"%d for d in ds))
print(" Gt   :", " ".join("%+7.1e"%tp[d] for d in ds))
print(" Gs   :", " ".join("%+7.1e"%sp[d] for d in ds))
with np.errstate(divide='ignore',invalid='ignore'): r=sp/tp
print(" Gs/Gt:", " ".join("%+7.3f"%r[d] for d in ds), "  (CFT: -2)")
PY
echo "===== DONE ====="
