#!/usr/bin/env bash
# run_glue_gevp_redo_claude.sh -- GEVP spectra for the massless redo ensembles (latest strategy:
# per-m GEVP, therm cut kmin=20). Uses the LOCAL analysis binary glue_gevp_analysis_claude.o
# (moved out of the obsolete src/both_3d 2026-08-13; takes the data dir as arg1).
# F   (l=1):  per_m=1, nops2=1, lsel="1"  -> 3 m-block grounds; jkdump for fit_perm.
# F^2 (0++):  per_m=0, nops2=2, l=0        -> 0++ = 2nd-lightest (index nops2-2); jkdump.
# Output in production/: gevp_{F,f2}_<tag>_claude.dat + _jk.dat + .log.  L4 skipped (too few cfgs).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
# ANA=../both_3d/glue_gevp_analysis_claude.o   # OLD: src/both_3d is OBSOLETE (NM 2026-08-13)
ANA=./glue_gevp_analysis_claude.o
if [ ! -x "$ANA" ]
then
  H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
  H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
  echo "### building $ANA ###"
  g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o "$ANA" \
    || { echo "### ANALYSIS BUILD FAILED ###"; exit 1; }
fi
AT=0.2
TCUT=14
BIN=0        # auto Nc/50
KMIN=20      # thermalization cut
RTOL=1e-8
TREBASE=1
NWORK=6

ENS=()
for nf in 2 4 6; do
  # L1, L2 (at0.2) + NEW L3 (gsq 1.5/3/4.5) + L4 (gsq 2/4/6) -- glue is 100% measured at every L
  for spec in "1 0.500000" "1 1.000000" "1 1.500000" \
              "2 1.000000" "2 2.000000" "2 3.000000" \
              "3 1.500000" "3 3.000000" "3 4.500000" \
              "4 2.000000" "4 4.000000" "4 6.000000"; do
    set -- $spec
    ENS+=("$nf $1 $2")   # nf L gsq
  done
done

ana_one () {
  local nf=$1 L=$2 g=$3
  local hb="hb1.000000"
  # L3 and L4 both use the resc-shifted hb tag
  { [ "$L" = "3" ] || [ "$L" = "4" ]; } && hb="hb0.400000-1.000000"
  local dd="data_Nf${nf}_gsq${g}at0.200000nu01.000000mRe0.000000mIm0.000000nt128L${L}_${hb}"
  [ -d "$dd" ] || { echo "skip (no dir) Nf${nf} g${g} L${L}"; return; }
  local tag="Nf${nf}_g${g}_L${L}"
  # F l=1 (per-m grounds) -- all L
  $ANA "$dd" 1 1 0 $AT gevp_F_${tag}_claude.dat 1 $TCUT $BIN $KMIN $RTOL glue_msm_shapes $TREBASE "" 0 1 0 "1" 1 0 0 gevp_F_${tag}_jk_claude.dat 1 > ana_F_${tag}_claude.log 2>&1
  # F l=2 (per-m grounds, H irrep; same F operator, lsel="2") -- ALL L now (L3 signal is good;
  # at L4 the low-stat ensembles may return nan, which the ratio scripts skip).
  $ANA "$dd" 1 1 0 $AT gevp_Fl2_${tag}_claude.dat 1 $TCUT $BIN $KMIN $RTOL glue_msm_shapes $TREBASE "" 0 1 0 "2" 1 0 0 gevp_Fl2_${tag}_jk_claude.dat 1 > ana_Fl2_${tag}_claude.log 2>&1
  # F^2 0++ (l=0) -- ALL L.  INCLUDES F^4 (NM 2026-08-18): the F^2-v2 basis has 14 shapes
  # (0-6 = p2 = F^2, 7-13 = p4 = F^4); the F^4 operators sharpen the 0++ variational estimate.
  # nops2=2 (NOT 4): keep only the 2 lightest GEVP modes = the near-zero VACUUM/constant mode
  # (vacsub=0) + the physical 0++.  The 0++ is then the NON-vacuum one = STATE 0 / cols (3,4),
  # ROBUST across all 36 ensembles by construction (no fixed high-index to misfire).
  # WHY NOT nops2=4: verify_f2v2_states_claude.py showed the fixed "0++ = state 2" assignment BREAKS
  # on the 5 noisiest ensembles (L4 g4/g6, L3 Nf4 g4.5): there state 2 catches a vacuum-like/spurious
  # mode and the 0++ falls to state 1.  nops2=4 err/ref 0.933 vs nops2=2 0.948 -- 1.5% traded for
  # robustness (glue_f2_v2_shapes_impl_plan_claude.md sec 6b: full_n2_v0 vs full_n4_v0).
  # OLD (7-shape F^2 only, no F^4): same nops2=2 / state 0 / cols (3,4), prefix glue_f2_shapes:
  # $ANA "$dd" 1 1 0 $AT gevp_f2_${tag}_claude.dat 2 $TCUT $BIN $KMIN $RTOL glue_f2_shapes $TREBASE "" 0 1 0 "0" 1 0 0 gevp_f2_${tag}_jk_claude.dat 0 > ana_f2_${tag}_claude.log 2>&1
  $ANA "$dd" 1 1 0 $AT gevp_f2_${tag}_claude.dat 2 $TCUT $BIN $KMIN $RTOL glue_f2_v2_shapes $TREBASE "" 0 1 0 "0" 1 0 0 gevp_f2_${tag}_jk_claude.dat 0 > ana_f2_${tag}_claude.log 2>&1
  echo "done $tag"
}
export -f ana_one
export ANA AT TCUT BIN KMIN RTOL TREBASE
printf '%s\n' "${ENS[@]}" | xargs -P "$NWORK" -I{} bash -c 'ana_one $1' _ {}
echo "ALL DONE"
