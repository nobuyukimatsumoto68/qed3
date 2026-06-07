#!/usr/bin/env bash
# jj_disc block-t validation: jj_disc_block_t == jj_disc (free field) via h5diff -d 1e-10 (NOT bit-exact:
# block K's Term B solve_shift_block differs ~1e-15). Validates BOTH the vector loops (massless) and the
# parity tilde loops (--mass-im 0.1). Both programs write the same disc.0.h5 per esnid, so run REF first,
# stash, run BLOCK_T, h5diff. Build included (fails fast on a compile error).
#
# Pick a FREE gpu:  GPU=0 bash tmp_claude3.sh   (default GPU 1).  jj_disc is cheap (apply-only, no
# connected source legs) so each run is a few min.
#
# Run from src/both_3d:   [GPU=0] bash tmp_claude3.sh
set -u
GPU=${GPU:-1}

LOG=disc_blockt_h5diff_claude.log
WARN=disc_blockt_build_warnings_claude.log
H5DIFF=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/bin/h5diff
REF=jj_disc_claude.o
BLK=jj_disc_block_t_claude.o

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "# jj_disc block-t h5diff (free field, GPU $GPU)  --  $(date)"

echo; echo "=== [0] build $REF + $BLK (warnings -> $WARN) ==="
if make -j4 "$REF" "$BLK" 2>"$WARN"; then
  echo "BUILD OK"
else
  echo "BUILD FAILED -- last 60 lines of $WARN:"; tail -60 "$WARN"; exit 1
fi

# one (label, esnid, extra-args) case per line
run_case(){
  local label="$1" esnid="$2"; shift 2
  local args="$*"
  local outdir="data_${esnid}/disc_nhits1"
  local h5="$outdir/disc.0.h5"
  local h5ref="$outdir/disc.0.ref.h5"
  echo; echo "=== CASE [$label]  args: $args ==="
  rm -f "$h5" "$h5ref"
  echo "-- REF $REF ($(date +%T)) --"
  if ! CUDA_VISIBLE_DEVICES=$GPU ./"$REF" $args --nhits 1; then echo "REF FAILED"; return 1; fi
  [ -f "$h5" ] || { echo "REF produced no $h5"; return 1; }
  mv "$h5" "$h5ref"
  echo "-- BLOCK_T $BLK ($(date +%T)) --"
  if ! CUDA_VISIBLE_DEVICES=$GPU ./"$BLK" $args --nhits 1; then echo "BLOCK_T FAILED"; return 1; fi
  [ -f "$h5" ] || { echo "BLOCK_T produced no $h5"; return 1; }
  echo "-- h5diff -d 1e-10 --"
  if "$H5DIFF" -d 1e-10 "$h5ref" "$h5"; then
    echo "[$label] H5DIFF within 1e-10 -- PASS"
  else
    echo "[$label] H5DIFF differences (>1e-10):"; "$H5DIFF" -v -d 1e-10 "$h5ref" "$h5" | head -60; echo "[$label] FAIL"
  fi
}

# vector (massless): exercises loops 1-2 (tp/sp/ylm).  parity (mass-im): adds loops 3-4 (tilde Jtil).
run_case "vector" "free_vmRe0.000000vmIm0.000000"
run_case "parity" "free_vmRe0.000000vmIm0.100000" --mass-im 0.1

echo; echo "=== done; log=$LOG  warnings=$WARN ==="
