#!/bin/bash -l
set -u

# Build + run the GEVP glueball-mass analysis on the h5 output of the massless L=1
# gsq8 ensembles. Writes one gnuplot-ready .dat per (operator, Nf):
#   col 1 = t   col 2 = ground-state mass   col 3 = error   [then per-state mean err]
# Plot e.g.:  plot 'gevp_f2_Nf2_claude.dat' u 1:2:3 w yerrorbars
#
# F^2  : N_FLOW=3 NCH=2 drop_l0=0  (full 96-op tower, l=0 kept)
# F_12 : N_FLOW=3 NCH=1 drop_l0=1  (45-op tower, l=0 dropped)
# Reads per-config F_corr.<k>.h5 (fast binary). Sequential over ensembles (I/O-bound).

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

AT=0.2
NOPS2=8      # GEVP states kept
TCUT=5       # effective mass at t = 1..tcut-1 (times a_t)
BINSIZE=0    # 0 => auto (~50 jackknife bins)
KMIN=100     # skip thermalization
RTOL=1e-8    # GEVP-metric pruning threshold (SVD prune of near-singular directions)

LOG=glue_gevp_run_claude.log

echo "================ BUILD analysis (g++ + h5) ================" | tee "$LOG"
date | tee -a "$LOG"
g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o glue_gevp_analysis_claude.o 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "BUILD FAILED rc=$rc" | tee -a "$LOG"; exit "$rc"; fi
echo "BUILD OK" | tee -a "$LOG"

dir_of () { echo "data_Nf${1}_gsq8.000000at0.200000nu01.000000nt128L1_${2}"; }

run_one () {
  local tag=$1 nf=$2 nch=$3 dropl0=$4
  local d
  d=$(dir_of "$nf" "$tag")
  local out="gevp_${tag}_Nf${nf}_claude.dat"
  echo "================ ${tag} Nf${nf} -> ${out} ================" | tee -a "$LOG"
  date | tee -a "$LOG"
  ./glue_gevp_analysis_claude.o "$d" 3 "$nch" "$dropl0" "$AT" "$out" "$NOPS2" "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" 2>&1 | tee -a "$LOG"
}

for NF in 2 4 6; do
  run_one f2  "$NF" 2 0
done
for NF in 2 4 6; do
  run_one msm "$NF" 1 1
done

echo "================ ALL DONE ================" | tee -a "$LOG"
date | tee -a "$LOG"
ls -la gevp_*_claude.dat 2>/dev/null | tee -a "$LOG"
