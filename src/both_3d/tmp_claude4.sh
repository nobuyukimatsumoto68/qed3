#!/usr/bin/env bash
# L1 (HMC force) exactness + gain: build + run test_grad_l1_claude.o on GPU 0 (jj free run is on GPU 1).
# Validates grad_deviceAsyncLaunch_l1 == grad_deviceAsyncLaunch per link (max|ref-L1| < 1e-10) and
# benchmarks the force link-sweep (grad vs grad_l1) -> speedup. Test calls BOTH methods explicitly, so
# build WITHOUT -DGRAD_L1.
#
# Run from src/both_3d:   bash tmp_claude4.sh
set -u
GPU=${GPU:-0}

LOG=l1_grad_test_claude.log
WARN=l1_build_warnings_claude.log
T=test_grad_l1_claude.o

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "# L1 grad exactness+gain test (GPU $GPU)  --  $(date)"

echo; echo "=== [0] build $T (warnings -> $WARN) ==="
if make -j4 "$T" 2>"$WARN"; then
  echo "BUILD OK"
else
  echo "BUILD FAILED -- last 60 lines of $WARN:"; tail -60 "$WARN"; exit 1
fi

echo; echo "=== [1] run ($(date +%T)) ==="
OMP_NUM_THREADS=4 CUDA_VISIBLE_DEVICES=$GPU ./"$T" || echo "test returned nonzero (FAIL)"

echo; echo "=== done; log=$LOG  warnings=$WARN ==="
