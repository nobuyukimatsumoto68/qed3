#!/usr/bin/env bash
# tmp_glue_therm_claude.sh -- build the flow-trajectory driver (L1/L2/L4 binaries) and run the
# FIRST TRAJECTORY CHECK on 6 representative LOCAL redo ensembles (Nf2, weak/strong gsq per L),
# sampled across each stream (early configs show thermalization, late ones the typical trajectory
# for setting c in t0^2 E_s(t0) = c).
# Run:  bash tmp_glue_therm_claude.sh          (output tee'd to tmp_glue_therm_claude.log)
# Writes: glue_therm_L{1,2,4}_claude.o (binaries), tmp_glue_therm_claude.log, and per ensemble
#   data_<basename>/therm_series_claude.dat  (text; one line per config: k + E_s at
#   t_fl = 0.0, 0.2, ..., 2.0). Re-run = resume (already-done k lines are skipped).
# No rm anywhere in this script. Configs are read strictly read-only.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=4

LOG=tmp_glue_therm_claude.log

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/'
SRC=glue_therm_series_claude.cu

{
  echo "##################################################################"
  echo "### build glue_therm_series (L1/L2/L4)  $(date)"
  echo "##################################################################"
} | tee "$LOG"

for L in 1 2 4
do
  BIN=glue_therm_L${L}_claude.o
  echo "### building $BIN ###" | tee -a "$LOG"
  $NVCC $NVCCBASE $INCLUDES -DN_REFINE_CLI=$L -o "$BIN" "$SRC" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]
  then
    echo "### ERROR: build failed for L=$L (exit $rc) -- aborting ###" | tee -a "$LOG"
    exit "$rc"
  fi
done

# ---- trajectory check: args = gsq Nf nu0 kmax_run kmin stride ens_dir tmax ----
# L1 stride 100 (~20 cfgs of 1999), L2 stride 50 (~20 of ~1000+), L4 stride 2 (~13 of ~25).
# tmax per ensemble by r = gsq/L (NM 2026-07-17): r=1.5 -> 1.2, r=1.0 -> 0.8, r=0.5 -> 0.4;
# eps = tmax/100, save every 5 steps (set in the driver).
run_check () {   # $1=L  $2=gsq  $3=stride  $4=tmax  $5=ens_dir
  local L=$1
  local gsq=$2
  local stride=$3
  local tmax=$4
  local ens=$5
  echo "##################################################################" | tee -a "$LOG"
  echo "### L$L gsq$gsq stride $stride tmax $tmax : $ens  $(date)" | tee -a "$LOG"
  echo "##################################################################" | tee -a "$LOG"
  ./glue_therm_L${L}_claude.o "$gsq" 2 1.0 100000 1 "$stride" "$ens"/ "$tmax" 2>&1 | tee -a "$LOG"
  local rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]
  then
    echo "### ERROR: run failed (exit $rc) -- aborting ###" | tee -a "$LOG"
    exit "$rc"
  fi
}

run_check 1 0.5 100 2.0 Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
run_check 1 1.0 100 2.0 Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
run_check 1 1.5 100 2.0 Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
run_check 2 1.0  50 2.0 Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
run_check 2 2.0  50 2.0 Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
run_check 2 3.0  50 2.0 Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
run_check 4 2.0   2 2.0 Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
run_check 4 4.0   2 2.0 Nf2_gsq4.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
run_check 4 6.0   2 2.0 Nf2_gsq6.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000

{
  echo "### sample output (first 5 lines of each .dat): ###"
  for d in data_Nf2_*_hb*
  do
    [ -f "$d"/therm_series_v6_claude.dat ] || continue
    echo "--- $d ---"
    head -10 "$d"/therm_series_v6_claude.dat
  done
  echo "### done  $(date) ###"
} | tee -a "$LOG"
