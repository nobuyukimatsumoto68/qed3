#!/usr/bin/env bash
# SMOKE TEST for the parity (m_P) refactor (audit I1/I2/I3).  Free U=1, m_P=0.1i, nhits=1, n-t0=2,
# spin x td=2 => 4 dilution patterns -- exactly the config that used to abort on pattern 2 (duplicate HDF5
# key).  parity is MASS-driven (no free_field guard), so this free run DOES exercise the (--)/disc blocks.
# Builds to a DISTINCT binary (jj_corr_dilute_mPtest.o) -> no clash with prod/cond/mF/mP binaries.
# ~15-25 min (1 free hit, 4 patterns, parity tilde solves).  Read back: mP_smoketest_claude.log
#
# PASS = (a) the run completes (no duplicate-key abort) AND (b) the per-hit .h5 has, at the h0/ LEVEL:
#   h0/tp/Vmm, h0/sp/Vmm, h0/disc/tp/Jtil, h0/disc/sp/Jtil, h0/s1/Vmm  (and NO h0/t0_0/tp/Vmm nesting).
#
# Usage: bash tmp_mP_smoketest_claude.sh [GPU]   (default 0)
# NOTE: if a COMPLETE corr.0.h0.h5 already exists in the output dir it will be SKIPPED (resume-gating) and
#   the test won't run -- in that case rm it yourself first (no rm in this script by house rule).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES="${1:-0}"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=mP_smoketest_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

OUT=data_free_vmRe0.000000vmIm0.100000/corr_dil_nt02_nhits1_s1td2/corr.0.h0.h5

echo "### compile -> jj_corr_dilute_mPtest.o  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute_mPtest.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE m_P=0.1i nhits=1 n-t0=2 spin+td2 (4 patterns; GPU $CUDA_VISIBLE_DEVICES)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute_mPtest.o --nhits 1 --n-t0 2 --spin-dilution --time-dilution 2 --mass-im 0.1 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### RUN FAILED ($st) -- I1 NOT fixed (or other error) ###" | tee -a "$LOG"; exit 1; fi

echo "### key-layout check on $OUT  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 - "$OUT" <<'PY' 2>&1 | tee -a "$LOG"
import sys, h5py
fn = sys.argv[1]
need_h0    = ['h0/tp/Vmm','h0/sp/Vmm','h0/disc/tp/Jtil','h0/disc/sp/Jtil',
              'h0/s1/Vmm','h0/s2/Vmm','h0/s3/Vmm','h0/tp/Vpp','h0/sp/Vpp']
must_absent= ['h0/t0_0/tp/Vmm','h0/t0_1/tp/Vmm','h0/t0_0/sp/Vmm']   # old per-origin nesting (I2)
with h5py.File(fn,'r') as f:
    ok = ('complete' in f)
    print('complete flag :', ok)
    for k in need_h0:
        p = (k+'/real') in f
        print('  present', k, ':', p); ok = ok and p
    for k in must_absent:
        a = (k+'/real') not in f
        print('  absent ', k, ':', a); ok = ok and a
    print('RESULT:', 'PASS -- I1/I2/I3 fix confirmed' if ok else 'FAIL -- see above')
    sys.exit(0 if ok else 2)
PY
echo "### smoke test done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
