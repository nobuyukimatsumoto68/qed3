#!/usr/bin/env bash
# C6f-c: validate jj_corr_block_t (t-blocked sink) == jj_corr_mrhs on FREE FIELD via h5diff -d 1e-10
# (NOT bit-exact: block K's Term B solve_shift_block differs ~1e-15 from op_K). Runs on the VACANT GPU 1
# (B1a is on GPU 0). Both programs write the same data_free_.../corr_nt01_nhits1/corr.0.h5, so run the
# REF first, stash it, run BLOCK_T, then h5diff. Build included so a compile error fails fast.
#
# Run from src/both_3d:   bash tmp_claude2.sh
set -u

LOG=c6fc_h5diff_claude.log
WARN=c6fc_build_warnings_claude.log
H5DIFF=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/bin/h5diff

REF=jj_corr_mrhs_claude.o
BLK=jj_corr_block_t_claude.o
NT0=1 ; NHITS=1
ARGS="--n-t0 $NT0 --nhits $NHITS"        # no --ens-dir => free field (U=1)

OUTDIR=data_free_vmRe0.000000vmIm0.000000/corr_nt0${NT0}_nhits${NHITS}
H5=$OUTDIR/corr.0.h5
H5REF=$OUTDIR/corr.0.ref.h5

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "# C6f-c jj_corr_block_t h5diff (free field, GPU 1)  --  $(date)"
echo "# args: $ARGS"

echo; echo "=== [0] build $REF + $BLK (warnings -> $WARN) ==="
if make -j4 "$REF" "$BLK" 2>"$WARN"; then
  echo "BUILD OK"
else
  echo "BUILD FAILED -- last 60 lines of $WARN:"; tail -60 "$WARN"; exit 1
fi

echo; echo "=== [1] clean stale outputs ==="
rm -f "$H5" "$H5REF"

echo; echo "=== [2] run REFERENCE  $REF  ($(date +%T)) ==="
if ! CUDA_VISIBLE_DEVICES=1 ./"$REF" $ARGS; then echo "REF run FAILED"; exit 1; fi
if [ ! -f "$H5" ]; then echo "REF produced no $H5"; exit 1; fi
mv "$H5" "$H5REF"
echo "# stashed reference -> $H5REF"

echo; echo "=== [3] run BLOCK_T  $BLK  ($(date +%T)) ==="
if ! CUDA_VISIBLE_DEVICES=1 ./"$BLK" $ARGS; then echo "BLOCK_T run FAILED"; exit 1; fi
if [ ! -f "$H5" ]; then echo "BLOCK_T produced no $H5"; exit 1; fi

echo; echo "=== [4] h5diff -d 1e-10  (expect: 0 differences within tol)  ($(date +%T)) ==="
if "$H5DIFF" -d 1e-10 "$H5REF" "$H5"; then
  echo "H5DIFF: within 1e-10 -- jj_corr_block_t == jj_corr_mrhs  PASS"
else
  rc=$?
  echo "H5DIFF: differences beyond 1e-10 (exit $rc) -- verbose detail:"
  "$H5DIFF" -v -d 1e-10 "$H5REF" "$H5" | head -100
  echo "H5DIFF FAIL"
fi

echo; echo "=== done; log=$LOG  warnings=$WARN ==="
