#!/usr/bin/env bash
# B1a (FIXED): thermalized-config h5diff -- jj_corr_mrhs == jj_corr BIT-IDENTICAL on ONE real config.
#
# BUG in the previous version: the ensemble has 3426 ckpoint_lat.* (dense; ckpoint_lat.448 exists), so
# --ninter 224 measured k=224,448,672,... (EVERY present multiple), ~46 min each => the REF loop ran for
# hours and never reached the h5diff. FIX: point --ens-dir at a temp dir holding a SINGLE symlink
# ckpoint_lat.224, so the loop measures k=224 then breaks at the missing k=448. One config, one hit.
#
# Both REF and MRHS use the same one-cfg dir => same esnid => one output path => h5diff. ~46 min PER
# program (~1.5 h total + build); the sink K applies (unchanged by mrhs) dominate each run.
#
# Run from src/both_3d:   bash tmp_claude.sh   (uses GPU 0 by default)
set -u

LOG=b1a_h5diff_claude.log
WARN=b1a_build_warnings_claude.log
H5DIFF=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/bin/h5diff

REALENS=Nf2_gsq8.000000at0.200000nu01.000000nt128L1_pole11
KCFG=224
ONECFG=b1a_onecfg_claude          # temp ens dir holding ONLY ckpoint_lat.224 (one symlink)
NINTER=$KCFG                      # loop: k=0 missing -> continue, k=224 -> measure, k=448 missing -> break
NHITS=1
NT0=1
GSQ=8 ; NF=2 ; NU0=1

REF=jj_corr_claude.o
MRHS=jj_corr_mrhs_claude.o

OUTDIR=data_${ONECFG}_vmRe0.000000vmIm0.000000/corr_nt0${NT0}_nhits${NHITS}
H5=$OUTDIR/corr.${KCFG}.h5
H5REF=$OUTDIR/corr.${KCFG}.ref.h5

ARGS="--gsq $GSQ --Nf $NF --nu0 $NU0 --ens-dir $ONECFG/ --ninter $NINTER --nhits $NHITS --n-t0 $NT0"

: > "$LOG"
exec > >(tee -a "$LOG") 2>&1
echo "# B1a thermalized h5diff (one-config)  --  $(date)"

echo; echo "=== [0] isolate a single config: $ONECFG/ckpoint_lat.$KCFG ==="
rm -rf "$ONECFG"; mkdir -p "$ONECFG"
ln -sfn "$PWD/$REALENS/ckpoint_lat.$KCFG" "$ONECFG/ckpoint_lat.$KCFG"
if [ ! -e "$ONECFG/ckpoint_lat.$KCFG" ]; then echo "symlink failed"; exit 1; fi
echo "# args: $ARGS"

echo; echo "=== [1] build ref + mrhs (warnings -> $WARN) ==="
if make -j4 "$REF" "$MRHS" 2>"$WARN"; then
  echo "BUILD OK"
else
  echo "BUILD FAILED -- last 50 lines of $WARN:"; tail -50 "$WARN"; exit 1
fi

echo; echo "=== [2] clean stale outputs for cfg $KCFG ==="
rm -f "$H5" "$H5REF"

echo; echo "=== [3] run REFERENCE  $REF  ($(date +%T)) ==="
if ! ./"$REF" $ARGS; then echo "REF run FAILED"; exit 1; fi
if [ ! -f "$H5" ]; then echo "REF produced no $H5"; exit 1; fi
mv "$H5" "$H5REF"
echo "# stashed reference -> $H5REF"

echo; echo "=== [4] run MRHS  $MRHS  ($(date +%T)) ==="
if ! ./"$MRHS" $ARGS; then echo "MRHS run FAILED"; exit 1; fi
if [ ! -f "$H5" ]; then echo "MRHS produced no $H5"; exit 1; fi

echo; echo "=== [5] h5diff  (expect: 0 differences, bit-identical)  ($(date +%T)) ==="
if "$H5DIFF" "$H5REF" "$H5"; then
  echo "H5DIFF: NO DIFFERENCES -- jj_corr_mrhs == jj_corr  BIT-IDENTICAL  PASS"
else
  rc=$?
  echo "H5DIFF: differences/error (exit $rc) -- verbose detail:"
  "$H5DIFF" -v "$H5REF" "$H5" | head -80
  echo "H5DIFF FAIL"
fi

echo; echo "=== done; log=$LOG  warnings=$WARN ==="
