#!/usr/bin/env bash
# run_glue_gevp_f2_v2_claude.sh -- Chunks 4 and 5 of glue_f2_v2_shapes_impl_plan_claude.md.
#
# PHASE A (VALIDATION GATE).  Analyse the new glue_f2_v2_shapes h5 restricted to shapes 0..6 -- the
# p=2 block, which IS the production F^2 basis -- and compare against the production glue_f2_shapes
# data analysed with byte-identical arguments.  Only the h5 prefix and orbits_arg differ, so the two
# must agree to round-off.  If they do not, the power axis broke something and PHASE B is skipped.
#
# PHASE B (full 14-operator basis).  orbits_arg empty = p=2 AND p=4 together, scanning nops2 and both
# vacsub settings, because the enlarged l=0 block carries more near-constant directions and the
# production bookkeeping (nops2=2, 0++ = state 0) will not survive unchanged.
#
# Analysis params are copied verbatim from run_glue_gevp_redo_claude.sh so the comparison is fair:
#   AT=0.2 TCUT=14 BIN=0(auto) KMIN=20 RTOL=1e-8 TREBASE=1, per_m=0, lsel="0", keepl0=1.
#
# Outputs (this dir): gevp_f2v2_<variant>_<tag>_claude.dat + _jk_claude.dat + ana_*.log
# Verdict table: glue_f2_v2_compare_claude.py (called at the end).
# ANALYSIS ONLY -- reads h5, writes .dat.  No rm anywhere.
#
#   bash run_glue_gevp_f2_v2_claude.sh 2>&1 | tee run_glue_gevp_f2_v2_claude.log
# Knobs: LEVEL, NF, GSQ, ENS, NOPS2_LIST, VACSUB_LIST
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}

# ANA=../both_3d/glue_gevp_analysis_claude.o   # OLD: src/both_3d is OBSOLETE (NM 2026-08-13)
ANA=./glue_gevp_analysis_claude.o
# build if absent (host-only g++ + Eigen + HighFive; recipe from the old both_3d glue_gevp_run script)
if [ ! -x "$ANA" ]
then
  H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
  H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
  echo "### building $ANA ###"
  g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o "$ANA" \
    || { echo "### ANALYSIS BUILD FAILED ###"; exit 1; }
fi
LEVEL="${LEVEL:-3}"
NF="${NF:-2}"
GSQ="${GSQ:-1.500000}"
ENS="${ENS:-Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000}"
NOPS2_LIST="${NOPS2_LIST:-2 3 4 5 6}"
VACSUB_LIST="${VACSUB_LIST:-0 1}"

# analysis constants -- identical to run_glue_gevp_redo_claude.sh
AT=0.2
TCUT=14
BIN=0
KMIN=20
RTOL=1e-8
TREBASE=1

PROD_PREFIX=glue_f2_shapes
V2_PREFIX=glue_f2_v2_shapes
PROD_ORBITS=""                 # production h5 has only the 7 shapes -> all of them
GATE_ORBITS="0,1,2,3,4,5,6"    # v2 h5: the p=2 block only
TAG="Nf${NF}_g${GSQ}_L${LEVEL}"
DD="data_${ENS}"

echo "================ F^2-v2 GEVP START $(date) ================"
echo "ens dir = $DD"
echo "tag     = $TAG"

# ---- preflight: the two prefixes must cover the SAME configs or the comparison is not fair ----
echo "---- preflight $(date) ----"
if [ ! -d "$DD" ]
then
  echo "ERROR: data dir not found: $DD"
  exit 1
fi
NPROD=$(ls "$DD"/${PROD_PREFIX}.*.h5 2>/dev/null | wc -l)
NV2=$(ls "$DD"/${V2_PREFIX}.*.h5 2>/dev/null | wc -l)
echo "production ${PROD_PREFIX} h5 : $NPROD"
echo "new        ${V2_PREFIX} h5 : $NV2"
if [ "$NV2" -eq 0 ]
then
  echo "ERROR: no ${V2_PREFIX} h5 -- run run_glue_f2_v2_shapes_claude.sh first"
  exit 1
fi
KPROD=$(ls "$DD"/${PROD_PREFIX}.*.h5 2>/dev/null | sed "s#.*/${PROD_PREFIX}\.##;s#\.h5##" | sort -n)
KV2=$(ls "$DD"/${V2_PREFIX}.*.h5 2>/dev/null | sed "s#.*/${V2_PREFIX}\.##;s#\.h5##" | sort -n)
if [ "$KPROD" = "$KV2" ]
then
  echo "config sets IDENTICAL -- gate comparison is apples-to-apples"
else
  echo "WARNING: config sets DIFFER (prod $NPROD vs v2 $NV2)."
  echo "         The gate compares two different samples; a mismatch would then NOT prove a bug."
  echo "         Finish the v2 sweep before trusting PHASE A."
fi

run_ana () {   # $1=prefix $2=orbits $3=nops2 $4=vacsub $5=variant-name
  local prefix="$1"
  local orbits="$2"
  local nops2="$3"
  local vacsub="$4"
  local name="$5"
  local out="gevp_f2v2_${name}_${TAG}_claude.dat"
  local jk="gevp_f2v2_${name}_${TAG}_jk_claude.dat"
  local log="ana_f2v2_${name}_${TAG}_claude.log"
  "$ANA" "$DD" 1 1 0 $AT "$out" "$nops2" $TCUT $BIN $KMIN $RTOL "$prefix" $TREBASE \
         "$orbits" "$vacsub" 1 0 "0" 1 0 0 "$jk" 0 > "$log" 2>&1
  local st=$?
  if [ "$st" -ne 0 ]
  then
    echo "  [FAIL status $st] $name  (see $log)"
    return 1
  fi
  if [ ! -s "$out" ]
  then
    echo "  [FAIL empty]  $name  (see $log)"
    return 1
  fi
  echo "  [ok] $name -> $out"
  return 0
}

# ================= PHASE A: validation gate =================
echo "---- PHASE A gate: p=2 sub-block vs production $(date) ----"
run_ana "$PROD_PREFIX" "$PROD_ORBITS" 2 0 "prod" || exit 1
run_ana "$V2_PREFIX"   "$GATE_ORBITS" 2 0 "gate" || exit 1

python3 glue_f2_v2_compare_claude.py gate \
  "gevp_f2v2_prod_${TAG}_claude.dat" "gevp_f2v2_gate_${TAG}_claude.dat"
GATE_ST=$?
if [ "$GATE_ST" -ne 0 ]
then
  echo "#### GATE FAILED -- the p=2 block does not reproduce production. PHASE B skipped. ####"
  echo "================ F^2-v2 GEVP STOPPED $(date) ================"
  exit 1
fi
echo "#### GATE PASSED ####"

# ================= PHASE B: full 14-operator basis =================
echo "---- PHASE B full basis (p=2 and p=4), nops2 x vacsub scan $(date) ----"
for vs in $VACSUB_LIST
do
  for n2 in $NOPS2_LIST
  do
    run_ana "$V2_PREFIX" "" "$n2" "$vs" "full_n${n2}_v${vs}" || true
  done
done

echo "---- verdict table $(date) ----"
python3 glue_f2_v2_compare_claude.py verdict "$TAG" "$NOPS2_LIST" "$VACSUB_LIST"
echo "================ F^2-v2 GEVP DONE $(date) ================"
