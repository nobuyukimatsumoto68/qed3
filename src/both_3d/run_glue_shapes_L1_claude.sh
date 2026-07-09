#!/bin/bash -l
set -u

# ============================================================================
# ONE-SHOT handoff: h5 generation (measurement) + GEVP analysis for the
# SHAPE-BASIS glueball operators on the massless L=1 gsq8 SEA ensembles (Nf2/4/6).
#
#   glue_f2_shapes_claude.cu    -> F^2 / 0++ (squared), ell 0-3   -> data_..._shapes_f2/
#   glue2_msm_shapes_claude.cu  -> linear F_12,          ell 1-3   -> data_..._shapes_msm/
#   glue_gevp_analysis_claude.cu (host g++ + Eigen + h5) reads the h5, GEVP -> .dat
#
# Basis = 5 icosahedral-orbit shapes {triangle, rect, twisted-rect, figure-8,
# twisted-figure-8}; single flow tmax=1.0/100; kmin=100 stride=1 (all configs).
# 6 single-threaded measurement processes (2 drivers x 3 Nf). Resume-safe (h5 complete flag).
# Analysis auto-detects nops from the h5, so it is L-agnostic. NO rm anywhere.
#
# *** BEFORE FIRST PRODUCTION RUN ***  the Nf2 data_..._shapes_{f2,msm}/ dirs hold 11
# STALE smoke h5 (old operator count). Resume checks only the complete flag, not nops,
# so a dir must be internally consistent -- rm the stale smoke dirs BY HAND first:
#   rm -r data_Nf2_gsq8.000000at0.200000nu01.000000nt128L1_shapes_f2
#   rm -r data_Nf2_gsq8.000000at0.200000nu01.000000nt128L1_shapes_msm
# (Nf4/Nf6 dirs are fresh.)
#
# Output: gevp_shapes_{f2,msm}_Nf{2,4,6}_claude.dat  (gnuplot: u 1:2:3 w yerrorbars)
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

# module load gcc/13.2.0
# module load cuda/12.6

export OMP_NUM_THREADS=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

GSQ=8.0
NU0=1.0
KMAXRUN=100000
KMIN=100
STRIDE=1
# analysis knobs
AT=0.2
NOPS2=8
TCUT=5
BINSIZE=0
RTOL=1e-8

LOG=glue_shapes_L1_claude.log

# ---------------------------------------------------------------------------
echo "================ BUILD (2 drivers + analysis) ================" | tee "$LOG"
date | tee -a "$LOG"
"$NVCC" glue2_msm_shapes_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_claude.o 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}; if [ "$rc" -ne 0 ]; then echo "linear BUILD FAILED rc=$rc" | tee -a "$LOG"; exit "$rc"; fi
"$NVCC" glue_f2_shapes_claude.cu   $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_claude.o   2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}; if [ "$rc" -ne 0 ]; then echo "F^2 BUILD FAILED rc=$rc" | tee -a "$LOG"; exit "$rc"; fi
g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o glue_gevp_analysis_claude.o 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}; if [ "$rc" -ne 0 ]; then echo "analysis BUILD FAILED rc=$rc" | tee -a "$LOG"; exit "$rc"; fi
echo "BUILD OK" | tee -a "$LOG"

# ---------------------------------------------------------------------------
echo "================ MEASURE (6 processes -> h5) ================" | tee -a "$LOG"
date | tee -a "$LOG"
for NF in 2 4 6; do
  ./glue_f2_shapes_claude.o   "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > glue_f2_shapes_Nf${NF}_run_claude.log  2>&1 &
  ./glue2_msm_shapes_claude.o "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > glue_msm_shapes_Nf${NF}_run_claude.log 2>&1 &
done
echo "launched 6 procs; waiting..." | tee -a "$LOG"
wait
echo "MEASURE DONE" | tee -a "$LOG"

# ---------------------------------------------------------------------------
echo "================ ANALYSIS (GEVP -> .dat) ================" | tee -a "$LOG"
date | tee -a "$LOG"
# tag  drop_l0 :  f2 (squared) keeps ell=0 ; msm (linear) drops ell=0
run_gevp () {
  local tag=$1 dropl0=$2 nf=$3
  local d="data_Nf${nf}_gsq8.000000at0.200000nu01.000000nt128L1_${tag}"
  local out="gevp_${tag}_Nf${nf}_claude.dat"
  echo "---- ${tag} Nf${nf} -> ${out} ----" | tee -a "$LOG"
  # N_FLOW/NCH (args 2,3) are ignored; nops is auto-detected from the h5.
  ./glue_gevp_analysis_claude.o "$d" 1 1 "$dropl0" "$AT" "$out" "$NOPS2" "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" 2>&1 | tee -a "$LOG"
}
for NF in 2 4 6; do run_gevp shapes_f2  0 "$NF"; done
for NF in 2 4 6; do run_gevp shapes_msm 1 "$NF"; done

echo "================ ALL DONE ================" | tee -a "$LOG"
date | tee -a "$LOG"
ls -la gevp_shapes_*_claude.dat 2>/dev/null | tee -a "$LOG"
