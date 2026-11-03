#!/usr/bin/env bash
# tmp_wilson_spec_claude.sh -- Wilson spectrum (paper Sec. IX.A.1): low eigenvalues of the bare Wilson D_W by
# shift-invert Arnoldi (eig_dw_arnoldi_prod_claude.cu, DW_ARNOLDI branch), m = 256 Krylov steps.
#   * all 31 Nf2 ensembles (12 massless at0.2, 4 massless at0.1, 15 massive), LAST config of each;
#   * shifts sigma = 0 (near the origin) and sigma = 1 (near the overlap wall M = 1);
#   * free field (cold gauge) per (L, a_t): at0.2 L1-L4, at0.1 L1-L2, both shifts.
# Output: final/analysis_wilson/spectrum/eig_dw_arnoldi_L<L>_<tag>_k<k>_sig<s>_nt128_claude.dat
# Plan: final/analysis_wilson/wilson_spectrum_impl_plan_claude.md
# GPU 0 only, one run at a time (each run is seconds to ~1 min).
# TOOLCHAIN (2026-09-30): /usr/local/cuda-12.9/bin/nvcc -ccbin g++-12; cuda-12.9 lib64 first.
#
# Run (detached):  nohup bash tmp_wilson_spec_claude.sh > /dev/null 2>&1 &
# Log:             tmp_wilson_spec_claude.log

set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=0

NM_ARNOLDI=256
NK=12
OUT_DIR=final/analysis_wilson/spectrum
mkdir -p "$OUT_DIR"

NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCBASE="-ccbin g++-12 -arch=sm_70 -O2 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/'
SRC=eig_dw_arnoldi_prod_claude.cu

LOG=tmp_wilson_spec_claude.log
echo "================ Wilson spectrum START $(date) ================" | tee -a "$LOG"

# ---- PHASE 1: build L1-L4 ----
echo "---- PHASE 1: build ----" | tee -a "$LOG"
for L in 1 2 3 4
do
  BIN=eig_dw_arnoldi_prod_L${L}_claude.o
  if [ -f "$BIN" ] && [ "$BIN" -nt "$SRC" ]
  then
    echo "### L${L}: $BIN up-to-date, skip ###" | tee -a "$LOG"
    continue
  fi
  echo "### compile L${L} -> $BIN [$(date '+%F %T')] ###" | tee -a "$LOG"
  $NVCC $NVCCBASE $INCLUDES -DLREF=$L -DNT_CLI=128 -o "$BIN" "$SRC" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]
  then
    echo "### ERROR: build failed for L=$L (rc=$rc) -- aborting ###" | tee -a "$LOG"
    exit "$rc"
  fi
done

# ---- PHASE 2: ensemble list (all Nf2) ----
ENS_LIST=$(ls -d Nf2_gsq*nt128L*_hb* 2>/dev/null)
NENS=$(printf '%s\n' "$ENS_LIST" | grep -c . || true)
echo "---- PHASE 2: $NENS Nf2 ensembles (expect 31) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_LIST" | tee -a "$LOG"
if [ "$NENS" -ne 31 ]
then
  echo "### ERROR: expected 31 Nf2 ensembles, found $NENS -- aborting ###" | tee -a "$LOG"
  exit 1
fi

NFAIL=0
NRUN=0

run_one () {
  # $1 = L, $2 = out tag, $3 = sigma ; CONFIG_LAT / AT_FREE set by the caller
  local L="$1"
  local tag="$2"
  local sig="$3"
  local dat="$OUT_DIR/eig_dw_arnoldi_L${L}${tag}_nt128_claude.dat"
  echo "[run $(date '+%F %T')] L=$L tag=$tag sigma=$sig" >> "$LOG"
  DW_ARNOLDI=1 DW_SHIFT="$sig" OUT_DIR="$OUT_DIR" OUT_TAG="$tag" \
    ./eig_dw_arnoldi_prod_L${L}_claude.o 0 0 -1 -1 12 "$NK" "$NM_ARNOLDI" >> "$LOG" 2>&1
  local rc=$?
  NRUN=$((NRUN+1))
  if [ "$rc" -ne 0 ] || [ ! -s "$dat" ]
  then
    echo "### FAILED: L=$L tag=$tag sigma=$sig (rc=$rc, dat present: $([ -s "$dat" ] && echo yes || echo no))" | tee -a "$LOG"
    NFAIL=$((NFAIL+1))
  fi
}

# ---- PHASE 3: interacting, last config, both shifts ----
echo "---- PHASE 3: interacting runs $(date) ----" | tee -a "$LOG"
for ens in $ENS_LIST
do
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  at=$(printf '%s' "$ens" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
  mre=$(printf '%s' "$ens" | grep -oE 'mRe[0-9.]+mIm' | sed 's/mRe//;s/mIm//')
  k=$(ls "$ens" | grep -E '^ckpoint_lat\.[0-9]+$' | sed 's/ckpoint_lat\.//' | sort -n | tail -1)
  for sig in 0 1
  do
    tag="_Nf2_gsq${gsq}_at${at}_mRe${mre}_k${k}_sig${sig}"
    CONFIG_LAT="$ens/ckpoint_lat.$k" run_one "$L" "$tag" "$sig"
  done
done

# ---- PHASE 4: free field per (L, a_t), both shifts ----
echo "---- PHASE 4: free-field runs $(date) ----" | tee -a "$LOG"
for spec in "1 0.2" "2 0.2" "3 0.2" "4 0.2" "1 0.1" "2 0.1"
do
  L=${spec% *}
  at=${spec#* }
  for sig in 0 1
  do
    tag="_free_at${at}_sig${sig}"
    AT_FREE="$at" run_one "$L" "$tag" "$sig"
  done
done

echo "---- SUMMARY: runs $NRUN, FAILED $NFAIL, dat files in $OUT_DIR: $(find "$OUT_DIR" -maxdepth 1 -name 'eig_dw_arnoldi_*_claude.dat' | wc -l) (expect 74) ----" | tee -a "$LOG"
echo "================ Wilson spectrum DONE $(date) ================" | tee -a "$LOG"
