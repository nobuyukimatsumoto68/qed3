#!/usr/bin/env bash
# run_sidewall_test_claude.sh -- chunk 2 of final/analysis_gluonic/jtop_psi_impl_plan_claude.md.
# Builds glue_sidewall_psi_claude.cu with -DSIDEWALL_TEST (L1) and runs the sanity checks on ONE config
# (Nf2 gsq1.0 L1 at0.2, k=20).  Test mode prints only (no h5 is written) and exits after that config.
# CPU-only, single thread.  NO rm, no kill.  All output -> sidewall_test_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}

LOG=sidewall_test_claude.log
: > "$LOG"

SRC=glue_sidewall_psi_claude.cu
BIN=glue_sidewall_psi_test_L1_claude.o
NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCFLAGS="-ccbin g++-12 -w -arch=sm_70 -O2 -std=c++17 -lcusolver -Xcompiler -fopenmp"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

echo "=== [1/2] build $BIN (-DSIDEWALL_TEST, L1) $(date) ===" | tee -a "$LOG"
"$NVCC" "$SRC" -DSIDEWALL_TEST $NVCCFLAGS $INCLUDES $LDFLAGS -o "$BIN" 2>&1 | tee -a "$LOG"
if [ "${PIPESTATUS[0]}" -ne 0 ] || [ ! -x "$BIN" ]
then
  echo "!! BUILD FAILED -- see $LOG" | tee -a "$LOG"
  exit 1
fi

ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
echo "" | tee -a "$LOG"
echo "=== [2/2] test run: $ENS, k=20 $(date) ===" | tee -a "$LOG"
# args: gsq Nf nu0 kmax_run kmin stride ens_dir
./"$BIN" 1.0 2 1.0 20 20 1 "$ENS" 2>&1 | tee -a "$LOG"
echo "=== done $(date) ===" | tee -a "$LOG"
