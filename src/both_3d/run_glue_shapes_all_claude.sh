#!/bin/bash -l
set -u

# ============================================================================
# COMPREHENSIVE handoff: shape-basis glueball measurement (both F and F^2) +
# GEVP analysis for ALL existing massless ensembles at L = 1, 2, 4 (all gsq,
# including the small-gsq scan). Glue h5 is written INTO THE SHARED data_ dir
# (same as the mesonic correlators), with distinct file prefixes:
#     glue_f2_shapes.<k>.h5    (F^2 / 0++, squared, ell 0-3)
#     glue_msm_shapes.<k>.h5   (linear F_12,        ell 1-3)
# Analysis -> gevp_{f2,msm}_Nf<nf>_gsq<gsq>_L<L>_claude.dat  (gnuplot u 1:2:3 w yerr).
#
# L is a COMPILE flag, so drivers are rebuilt once per L. Basis = 5 icosahedral-orbit
# shapes {triangle, rect, twisted-rect, figure-8, twisted-figure-8}; single flow 1.0/100.
# Resume-safe (h5 complete flag). NO rm anywhere.
#
# *** SCALE WARNING ***  nops grows with L (L1=80, L2=256, L4=784), so F_corr files are
# ~6.5 MB (L1) / 67 MB (L2) / 630 MB (L4) per config. Total output is HUNDREDS OF GB and
# the run takes HOURS. Adjust LLIST / STRIDE / KMIN below to taste before launching.
#
# Knobs (env-overridable):  LLIST="1 2 4"  NPROC=6  KMIN=20  STRIDE=1
#   KMIN=20 is a modest thermalization cut safe for the smallest ensemble (68 cfg). Raise it
#   for the long ensembles if you want a heavier cut (it is ensemble-independent here).
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

# module load gcc/13.2.0
# module load cuda/12.6
export OMP_NUM_THREADS=1

LLIST="${LLIST:-1 2 4}"
NPROC="${NPROC:-8}"   # one worker per ensemble; 8-way ensemble-level parallelism (OMP_NUM_THREADS=1 each)
KMIN="${KMIN:-20}"
STRIDE="${STRIDE:-1}"
NU0=1.0
KMAXRUN=100000
# analysis knobs
AT=0.2
NOPS2=8
TCUT=5
BINSIZE=0
RTOL=1e-8

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=glue_shapes_all_claude.log
echo "================ START $(date) ================" | tee "$LOG"

# build the (L-agnostic) analysis once
g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o glue_gevp_analysis_claude.o 2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "analysis BUILD FAILED" | tee -a "$LOG"; exit 1; }

# wait until fewer than NPROC background jobs are running
pool_wait () { while [ "$(jobs -rp | wc -l)" -ge "$NPROC" ]; do wait -n; done; }

for L in $LLIST; do
  echo "================ L=$L : BUILD DRIVERS ================" | tee -a "$LOG"
  DEFS=""; [ "$L" != "1" ] && DEFS="-DN_REFINE_CLI=$L"
  "$NVCC" glue2_msm_shapes_claude.cu $NVCCFLAGS $DEFS $INCLUDES $LDFLAGS -o glue2_msm_shapes_claude.o 2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "L=$L linear BUILD FAILED" | tee -a "$LOG"; exit 1; }
  "$NVCC" glue_f2_shapes_claude.cu   $NVCCFLAGS $DEFS $INCLUDES $LDFLAGS -o glue_f2_shapes_claude.o   2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "L=$L F^2 BUILD FAILED" | tee -a "$LOG"; exit 1; }

  # enumerate massless config dirs for this L
  ens=()
  for d in Nf*_gsq*at0.200000nu01.000000nt128L${L}; do
    [ -d "$d" ] || continue
    case "$d" in *mRe*|*mIm*) continue;; esac
    [ "$(ls "$d"/ckpoint_lat.* 2>/dev/null | wc -l)" -gt 0 ] && ens+=("$d")
  done
  echo "L=$L : ${#ens[@]} massless ensembles: ${ens[*]}" | tee -a "$LOG"

  echo "================ L=$L : MEASURE (one ensemble per worker, pool=$NPROC) ================" | tee -a "$LOG"
  for d in "${ens[@]}"; do
    nf=$(echo "$d"  | sed 's/^Nf\([0-9]*\)_.*/\1/')
    gsq=$(echo "$d" | sed 's/.*gsq\([0-9.]*\)at.*/\1/')
    pool_wait
    # one worker = one ensemble: both ops (F^2 then linear F_12) run sequentially inside the subshell
    (
      ./glue_f2_shapes_claude.o   "$gsq" "$nf" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > "glue_f2_${d}_run_claude.log"  2>&1
      ./glue2_msm_shapes_claude.o "$gsq" "$nf" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > "glue_msm_${d}_run_claude.log" 2>&1
    ) &
  done
  wait
  echo "L=$L MEASURE DONE $(date)" | tee -a "$LOG"

  echo "================ L=$L : ANALYSIS (GEVP -> .dat) ================" | tee -a "$LOG"
  for d in "${ens[@]}"; do
    nf=$(echo "$d"  | sed 's/^Nf\([0-9]*\)_.*/\1/')
    gsq=$(echo "$d" | sed 's/.*gsq\([0-9.]*\)at.*/\1/')
    # args 2,3 (N_FLOW,NCH) ignored (nops auto-detected); arg12 = h5 file prefix
    ./glue_gevp_analysis_claude.o "data_$d" 1 1 0 "$AT" "gevp_f2_Nf${nf}_gsq${gsq}_L${L}_claude.dat"  "$NOPS2" "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" glue_f2_shapes  2>&1 | tee -a "$LOG"
    # drop_l0=0 : keep ell=0 (monopole/topological) as an extra vacuum-subtracted variational op
    ./glue_gevp_analysis_claude.o "data_$d" 1 1 0 "$AT" "gevp_msm_Nf${nf}_gsq${gsq}_L${L}_claude.dat" "$NOPS2" "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" glue_msm_shapes 2>&1 | tee -a "$LOG"
  done
done

echo "================ ALL DONE $(date) ================" | tee -a "$LOG"
ls -la gevp_*_claude.dat 2>/dev/null | tee -a "$LOG"
