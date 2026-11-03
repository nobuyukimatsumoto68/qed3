#!/usr/bin/env bash
# build_gevp_uniform_claude.sh -- build the GEVP analysis binary with UNIFORM metric pruning (NM 2026-10-12:
# kept dimension decided once on the full sample, replayed in every jackknife sample) and A/B-check it against
# the old per-sample-pruning binary on one ensemble.  New binary: glue_gevp_analysis_uni_claude.o
# (the old glue_gevp_analysis_t0_claude.o is NOT touched).  Host-only g++ build, single thread.  NO rm.
# All output -> build_gevp_uniform_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
LOG=build_gevp_uniform_claude.log
: > "$LOG"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
NEW=./glue_gevp_analysis_uni_claude.o
OLD=./glue_gevp_analysis_t0_claude.o

echo "=== [1/2] build $NEW $(date) ===" | tee -a "$LOG"
g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include -I/opt/eigen-3.4.0/ $H5I glue_gevp_analysis_claude.cu $H5L -o "$NEW" 2>&1 | tee -a "$LOG"
if [ "${PIPESTATUS[0]}" -ne 0 ] || [ ! -x "$NEW" ]
then
  echo "!! BUILD FAILED -- see $LOG" | tee -a "$LOG"
  exit 1
fi

echo "" | tee -a "$LOG"
echo "=== [2/2] A/B check: old vs new on Nf2 g1.0 L1 at0.2, F l=1, binsize 80 $(date) ===" | tee -a "$LOG"
DD=data_Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
OUT=build_gevp_uniform_check_claude
mkdir -p "$OUT"
for B in OLD NEW
do
  if [ "$B" = "OLD" ]
  then
    ANA=$OLD
  else
    ANA=$NEW
  fi
  "$ANA" "$DD" 1 1 0 1.0 "$OUT/F_${B}_claude.dat" 1 14 80 20 1e-8 glue_msm_shapes 2 "" 0 1 0 "1" 1 0 0 \
    "$OUT/F_${B}_jk_claude.dat" 1 1 0 "0" > "$OUT/F_${B}_claude.log" 2>&1
  echo "$B exit=$?" | tee -a "$LOG"
done
grep -E "uniform pruning" "$OUT/F_NEW_claude.log" | tee -a "$LOG"
python3 - <<'PY' 2>&1 | tee -a "$LOG"
import numpy as np
a = np.loadtxt("build_gevp_uniform_check_claude/F_OLD_jk_claude.dat")
b = np.loadtxt("build_gevp_uniform_check_claude/F_NEW_jk_claude.dat")
print("jk shapes", a.shape, b.shape, " max|new-old| =", np.nanmax(np.abs(a - b)), " nan old/new:", int(np.isnan(a).sum()), int(np.isnan(b).sum()))
PY
echo "=== done $(date) ===" | tee -a "$LOG"
