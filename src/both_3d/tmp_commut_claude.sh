#!/usr/bin/env bash
# tmp_commut_claude.sh -- build + validate the commutator-tp current (jj_commut_tp_claude.cc) on the
# EXISTING free L1 propagator (data_free_.../prop_deter_L1/Dinv.0.h5 already carries Dov + Dm_inv).
# Pure host g++ (no CUDA, no GPU).  Computes the tp correlator via the [D_ov,Theta] collapse and prints
# Gt_avg(dt) + the disc (= P.D_ov - I diagonal) sanity check.
# Run:  bash tmp_commut_claude.sh 2>&1 | tee commut_tp_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load gcc/13.2.0 2>/dev/null

INCLUDES='-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5 -lm'

PROP=data_free_vmRe0.000000vmIm0.000000/prop_deter_L1/Dinv.0.h5
[ -f "$PROP" ] || { echo "missing $PROP"; exit 1; }

echo "===== compile jj_commut_tp.o (host g++) ====="
g++ -std=c++17 -O2 $INCLUDES jj_commut_tp_claude.cc -o jj_commut_tp.o $LDFLAGS \
  || { echo "BUILD FAILED"; exit 1; }

echo "===== run commutator-tp (free, L1, massless) ====="
rm -f data_free_vmRe0.000000vmIm0.000000/corr_deter_commut_L1/corr.0.h5*
./jj_commut_tp.o --mass-re 0.0 --mass-im 0.0 --L 1 --n-t0 2 \
  || { echo "RUN FAILED"; exit 1; }

echo "===== Gt_avg(dt) full profile (Python) ====="
OMP_NUM_THREADS=1 python3 - << 'PY'
import h5py, numpy as np
f=h5py.File('data_free_vmRe0.000000vmIm0.000000/corr_deter_commut_L1/corr.0.h5','r')
g=f['h0/tp/Gt_avg/real'][:]+1j*f['h0/tp/Gt_avg/imag'][:]
Nt=int(f['Nt'][0])
print(" dt   :", " ".join("%6d"%d for d in [1,2,3,4,6,8,12,16,24,32,48,64]))
print(" Gt(dt):", " ".join("%+.2e"%g[d].real for d in [1,2,3,4,6,8,12,16,24,32,48,64]))
print(" large-dt plateau (Q^2 pedestal?): Gt[48..80] ~", np.mean(g[48:81].real))
print(" Gt[1]/Gt[2]=%.3f  Gt[2]/Gt[4]=%.3f (decay shape)"%(g[1].real/g[2].real, g[2].real/g[4].real))
PY
echo "===== DONE ====="
