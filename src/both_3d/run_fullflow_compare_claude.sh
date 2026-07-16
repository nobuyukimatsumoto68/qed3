#!/bin/bash -l
set -u

# ============================================================================
# FULL-3D-FLOW comparison run (temporary, for A/B vs the spatial-only flow).
# Builds the F and F^2 shape drivers with -DFLOW_FULL so the Wilson gradient flow is driven by the
# FULL 3D gauge-action gradient (get_force = spatial + temporal plaquettes) instead of the spatial-only
# gradient (get_spatial, the "on a timeslice" restriction). Outputs go to a DISTINCT prefix
#   glue_msm_shapes_fullflow.*.h5   /   glue_f2_shapes_fullflow.*.h5
# in the SAME data_<cfgdir>/ dirs, so the existing spatial-flow production h5 are untouched.
#
# Comparison ensemble: Nf2 gsq8 at L1 and L2 (L2 F^2 is the marginal-signal case of interest).
# Runtime args per driver: gsq Nf nu0 kmax kmin stride.
# 4 single-threaded workers (F/F^2 x L1/L2) => 4 cores. FLOW_TMAX/NSTEP as in production (2.0 / 200).
#
# NOTE: this measures the SAME configs as production, so error bars are directly comparable. To do a
# quick look first, raise STRIDE (e.g. 4) for fewer configs. No rm is placed here; the _fullflow prefix
# is fresh, and the "complete" flag makes re-runs resume-safe.
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=1

NU0=1.0
GSQ=8.000000
NF=2
KMAX=100000
KMIN=1
STRIDE=1        # set to 4 for a fast preview (~1/4 configs)

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=fullflow_compare_claude.log
echo "================ FULL-FLOW compare build+run START $(date) ================" | tee "$LOG"

# ---- build 4 binaries with -DFLOW_FULL (L via N_REFINE_CLI) --------------------
echo "---- build $(date) ----" | tee -a "$LOG"
for L in 1 2; do
  if [ "$L" = "1" ]; then DEF="-DFLOW_FULL"; else DEF="-DFLOW_FULL -DN_REFINE_CLI=$L"; fi
  "$NVCC" glue2_msm_shapes_claude.cu $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_fullflow_L${L}_claude.o 2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F  L$L FULLFLOW BUILD FAILED" | tee -a "$LOG"; exit 1; }
  "$NVCC" glue_f2_shapes_claude.cu   $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_fullflow_L${L}_claude.o   2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F^2 L$L FULLFLOW BUILD FAILED" | tee -a "$LOG"; exit 1; }
done
echo "---- build OK $(date) ----" | tee -a "$LOG"

# ---- run: 4 workers (F/F^2 x L1/L2) on Nf2 gsq8 -------------------------------
run_one() {
  local chan=$1
  local L=$2
  local bin elog
  if [ "$chan" = "F" ]; then bin=./glue2_msm_shapes_fullflow_L${L}_claude.o; else bin=./glue_f2_shapes_fullflow_L${L}_claude.o; fi
  elog="fullflow_${chan}_L${L}_claude.log"
  echo "[start $(date '+%F %T')] $chan L$L" >> "$LOG"
  "$bin" "$GSQ" "$NF" "$NU0" "$KMAX" "$KMIN" "$STRIDE" > "$elog" 2>&1
  echo "[done  $(date '+%F %T')] $chan L$L  (log: $elog)" >> "$LOG"
}
export -f run_one
export GSQ NF NU0 KMAX KMIN STRIDE LOG

echo "---- launch 4-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "F 1" "F2 1" "F 2" "F2 2" | xargs -P 4 -I{} bash -c 'run_one $1 $2' _ {}

echo "================ FULL-FLOW compare DONE $(date) ================" | tee -a "$LOG"
