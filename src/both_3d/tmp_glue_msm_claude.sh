#!/bin/bash -l
set -u

# Build (and optionally run) glue2_msm_claude.cu -- items 2 (norm fix) + 7 (multi-flow).
# Mirrors the Makefile NVCC flags + include paths (qfe_mod/include for Eigen,
# includes/ for project headers). Adjust paths if your working build differs.
# Run from src/both_3d/. Output object: glue2_msm_claude.o
#
# Usage:
#   bash tmp_glue_msm_claude.sh            # build only
#   RUN=1 bash tmp_glue_msm_claude.sh      # build, then run with default args (gsq=2 Nf=2 nu0=1)
#   RUN=1 GSQ=2.0 NF=2 NU0=1.0 bash tmp_glue_msm_claude.sh

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

LOG=glue2_msm_build_claude.log
NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/"

# module load gcc/13.2.0   # uncomment on the cluster
# module load cuda/12.6

echo "================ BUILD glue2_msm_claude.cu ================" | tee "$LOG"
date | tee -a "$LOG"
"$NVCC" glue2_msm_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_claude.o 2>&1 | tee -a "$LOG"
build_rc=${PIPESTATUS[0]}
if [ "$build_rc" -ne 0 ]; then
  echo "BUILD FAILED (rc=$build_rc) -- see $LOG" | tee -a "$LOG"
  exit "$build_rc"
fi
echo "BUILD OK -> glue2_msm_claude.o" | tee -a "$LOG"

if [ "${RUN:-0}" = "1" ]; then
  GSQ=${GSQ:-2.0}
  NF=${NF:-2}
  NU0=${NU0:-1.0}
  echo "================ RUN glue2_msm_claude.o $GSQ $NF $NU0 ================" | tee -a "$LOG"
  echo "(writes correlators to data_*_msm/ ; reads ckpoints from the un-suffixed dir)" | tee -a "$LOG"
  date | tee -a "$LOG"
  ./glue2_msm_claude.o "$GSQ" "$NF" "$NU0" 2>&1 | tee -a "$LOG"
  run_rc=${PIPESTATUS[0]}
  echo "RUN rc=$run_rc" | tee -a "$LOG"
fi
