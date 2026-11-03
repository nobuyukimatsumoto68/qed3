#!/usr/bin/env bash
# tmp_flow_move_claude.sh -- move the Euler-draft flow_weighted_claude.h5 (401-point grid, tmax 4.0) out of the
# way so the RK3 run can start, then launch tmp_flow_weighted_claude.sh. Moves only (mv -n, never overwrites).
#   1. rename the misnamed backup old_flow_weighted_euler_claude/_flow_weighted_claude.h5 (= L1 g0.5 Euler, 3999 cfg)
#   2. move the 11 remaining Euler-draft files (L1 g0.5 is NOT in the list: its current file is the new RK3 one)
#   3. launch the flow run (detached)
# Run:  bash tmp_flow_move_claude.sh      Log: tmp_flow_move_claude.log

set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
LOG=tmp_flow_move_claude.log
BK=old_flow_weighted_euler_claude
echo "==== move START $(date) ====" | tee -a "$LOG"

# ---- 1. misnamed backup ----
SRC1="$BK/_flow_weighted_claude.h5"
DST1="$BK/Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000_flow_weighted_claude.h5"
if [ -f "$SRC1" ]
then
  mv -n "$SRC1" "$DST1"
  if [ -f "$SRC1" ]
  then
    echo "### ERROR: could not rename $SRC1 -- aborting ###" | tee -a "$LOG"
    exit 1
  fi
  echo "renamed $SRC1 -> $DST1" | tee -a "$LOG"
fi

# ---- 2. the 11 remaining Euler-draft files ----
ENS11="
Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000
Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000
Nf2_gsq4.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000
Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
Nf2_gsq4.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
Nf2_gsq6.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
"
for ens in $ENS11
do
  src="data_${ens}/flow_weighted_claude.h5"
  dst="$BK/${ens}_flow_weighted_claude.h5"
  if [ ! -f "$src" ]
  then
    echo "skip (no file): $src" | tee -a "$LOG"
    continue
  fi
  if [ -e "$dst" ]
  then
    echo "### ERROR: $dst already exists -- not overwriting; aborting ###" | tee -a "$LOG"
    exit 1
  fi
  mv -n "$src" "$dst"
  if [ -f "$src" ]
  then
    echo "### ERROR: mv failed for $src -- aborting ###" | tee -a "$LOG"
    exit 1
  fi
  echo "moved $src -> $dst" | tee -a "$LOG"
done
echo "backup dir now holds $(find "$BK" -maxdepth 1 -name '*.h5' | wc -l) files (expect 12)" | tee -a "$LOG"

# ---- 3. launch the RK3 flow run ----
nohup bash tmp_flow_weighted_claude.sh > /dev/null 2>&1 &
echo "launched tmp_flow_weighted_claude.sh (pid $!) $(date)" | tee -a "$LOG"
echo "==== move DONE ====" | tee -a "$LOG"
