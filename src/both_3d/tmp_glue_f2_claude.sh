#!/bin/bash -l
set -u

# Build + smoke-run glue_f2_claude.cu (F_{\mu\nu}^2 action-density / 0^{++} glueball).
# Mirrors the Makefile NVCC flags + include paths. Run from src/both_3d/.
# Refinement L is a COMPILE flag (-DN_REFINE_CLI); default L=1.
#
# Usage:
#   bash tmp_glue_f2_claude.sh                    # build (L=1) + smoke run (Nf2 gsq8, 5 configs)
#   L=2 bash tmp_glue_f2_claude.sh                # build L=2 + smoke run
#   RUN=0 bash tmp_glue_f2_claude.sh              # build only
#   L=1 GSQ=8.0 NF=2 NU0=1.0 NCONF=5 bash tmp_glue_f2_claude.sh

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

L=${L:-1}
GSQ=${GSQ:-8.0}
NF=${NF:-2}
NU0=${NU0:-1.0}
NCONF=${NCONF:-5}
RUN=${RUN:-1}

LOG=glue_f2_smoke_claude.log
BIN=glue_f2_claude.o
NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/"
DEFS=""
if [ "$L" != "1" ]; then
  DEFS="-DN_REFINE_CLI=$L"
fi

# module load gcc/13.2.0   # uncomment on the cluster
# module load cuda/12.6

echo "================ BUILD glue_f2_claude.cu (L=$L) ================" | tee "$LOG"
date | tee -a "$LOG"
"$NVCC" glue_f2_claude.cu $NVCCFLAGS $DEFS $INCLUDES $LDFLAGS -o "$BIN" 2>&1 | tee -a "$LOG"
build_rc=${PIPESTATUS[0]}
if [ "$build_rc" -ne 0 ]; then
  echo "BUILD FAILED (rc=$build_rc) -- see $LOG" | tee -a "$LOG"
  exit "$build_rc"
fi
echo "BUILD OK -> $BIN" | tee -a "$LOG"

if [ "$RUN" = "1" ]; then
  # match C++ std::to_string formatting (6 decimals) for the dir names
  GSQF=$(printf '%.6f' "$GSQ")
  NU0F=$(printf '%.6f' "$NU0")
  DIR3="Nf${NF}_gsq${GSQF}at0.200000nu0${NU0F}nt128L${L}"
  DIR4="data_${DIR3}_f2"
  echo "================ SMOKE RUN $BIN  gsq=$GSQ Nf=$NF nu0=$NU0  ($NCONF configs) ================" | tee -a "$LOG"
  echo "(reads ckpoints from ${DIR3}/ ; writes ${DIR4}/)" | tee -a "$LOG"
  date | tee -a "$LOG"
  ./"$BIN" "$GSQ" "$NF" "$NU0" "$NCONF" 2>&1 | tee -a "$LOG"
  run_rc=${PIPESTATUS[0]}
  echo "RUN rc=$run_rc" | tee -a "$LOG"
  echo "---- output dir listing ----" | tee -a "$LOG"
  ls -la "${DIR4}/" 2>&1 | tee -a "$LOG"
fi
