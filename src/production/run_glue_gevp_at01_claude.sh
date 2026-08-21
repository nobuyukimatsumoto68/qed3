#!/usr/bin/env bash
# run_glue_gevp_at01_claude.sh -- GEVP F (l=1) + Fl2 for the at=0.1 PAIRED ensembles (L1 g1.0, L2 g2.0,
# Nf2/4/6), tagged _at0.1 (AT=0.1) so the at=0.2 gevp files are NOT clobbered.  For the a_t^2 figure.
# Same analysis binary / parameters as run_glue_gevp_redo_claude.sh; F l=1 uses glue_msm_shapes.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
ANA=./glue_gevp_analysis_claude.o
[ -x "$ANA" ] || { echo "### missing $ANA -- build via run_glue_gevp_redo first ###"; exit 1; }
AT=0.1          # half-a_t ensembles
TCUT=14
BIN=0
KMIN=20
RTOL=1e-8
TREBASE=1
NWORK=4

# paired at0.1 ensembles: (L gsq) = (1 1.0), (2 2.0)
ENS=()
for nf in 2 4 6; do
  for spec in "1 1.000000" "2 2.000000"; do
    set -- $spec
    ENS+=("$nf $1 $2")
  done
done

ana_one () {
  local nf=$1 L=$2 g=$3
  local dd="data_Nf${nf}_gsq${g}at0.100000nu01.000000mRe0.000000mIm0.000000nt128L${L}_hb1.000000"
  [ -d "$dd" ] || { echo "skip (no dir) Nf${nf} g${g} L${L} at0.1"; return; }
  local tag="Nf${nf}_g${g}_L${L}_at0.1"
  $ANA "$dd" 1 1 0 $AT gevp_F_${tag}_claude.dat 1 $TCUT $BIN $KMIN $RTOL glue_msm_shapes $TREBASE "" 0 1 0 "1" 1 0 0 gevp_F_${tag}_jk_claude.dat 1 > ana_F_${tag}_claude.log 2>&1
  $ANA "$dd" 1 1 0 $AT gevp_Fl2_${tag}_claude.dat 1 $TCUT $BIN $KMIN $RTOL glue_msm_shapes $TREBASE "" 0 1 0 "2" 1 0 0 gevp_Fl2_${tag}_jk_claude.dat 1 > ana_Fl2_${tag}_claude.log 2>&1
  echo "done $tag (dir $dd)"
}
export -f ana_one
export ANA AT TCUT BIN KMIN RTOL TREBASE
printf '%s\n' "${ENS[@]}" | xargs -P "$NWORK" -I{} bash -c 'ana_one $1' _ {}
echo "ALL DONE at0.1 GEVP"
