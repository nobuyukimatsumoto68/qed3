#!/usr/bin/env bash
# tmp_flow_topup_claude.sh -- Wilson-flow TOP-UP / REDO (NM 2026-10-24), 30 ensembles:
#   (a) 18 STALE massless at0.2 L3/L4 ensembles (flow h5 dated 2026-08-08; configs grew since:
#       L3 -> 999, L4 -> 799). Resume flows ONLY the missing configs.
#   (b) 12 massless at0.1 ensembles (L1 g0.5/g1.0/g1.5, L2 g2.0; Nf2/4/6), all from scratch:
#       - 6 (L1 g1.0, L2 g2.0) had flow h5 made with the OLD at=0.2 hardcode (wrong beta_s).
#         PHASE 0 MOVES them to old_at0.1_flow_claude/ (kept for comparison, nothing deleted),
#         so the driver starts a fresh file with the correct a_t.
#       - 6 (L1 g0.5, g1.5) never had a flow h5.
# Massive ensembles are NOT included: their flow h5 already cover every config.
# Protocol v3 unchanged: tmax = 20 gsq/L, 400 steps, save every 2 (201 points).
# Driver: glue_therm_flow_claude.cu (a_t now read from the ens dir; at0.2 -> 0.2 = unchanged,
# at0.1 -> 0.1). CPU-only, single-threaded workers, one ensemble per worker.
# TOOLCHAIN (2026-09-30): /usr/local/cuda-12.9/bin/nvcc -ccbin g++-12; cuda-12.9 lib64 first.
# New binary names (_cu129) so the old 08-08 binaries are left untouched.
#
# Run (detached):  nohup bash tmp_flow_topup_claude.sh > /dev/null 2>&1 &
# Log:             tmp_flow_topup_claude.log   (per-ensemble: therm_flow_<tag>_topup_claude.log)
# Workers:         NWORK=12 by default (override: NWORK=8 nohup bash tmp_flow_topup_claude.sh ...)

set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

NWORK=${NWORK:-12}
KMAX=100000
KMIN=1
STRIDE=1

NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCBASE="-ccbin g++-12 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lm'
SRC=glue_therm_flow_claude.cu

PAT_STALE='Nf*_gsq*at0.200000nu01.000000mRe0.000000mIm0.000000nt128L[34]_hb0.400000-1.000000'
PAT_AT01='Nf*_gsq*at0.100000nu01.000000mRe0.000000mIm0.000000nt128L[12]_hb1.000000'
BACKUP=old_at0.1_flow_claude

LOG=tmp_flow_topup_claude.log
echo "================ flow top-up START $(date) ================" | tee -a "$LOG"

# ---- PHASE 0: move the 6 wrong-a_t at0.1 flow h5 aside (no delete; abort if a backup exists) ----
echo "---- PHASE 0: move old at0.1 flow h5 -> $BACKUP/ ----" | tee -a "$LOG"
mkdir -p "$BACKUP"
for d in $PAT_AT01
do
  [ -d "$d" ] || continue
  src="data_$d/therm_flow_claude.h5"
  [ -f "$src" ] || continue
  dst="$BACKUP/${d}_therm_flow_claude.h5"
  if [ -e "$dst" ]
  then
    echo "### ERROR: $dst already exists -- not overwriting a backup; aborting ###" | tee -a "$LOG"
    exit 1
  fi
  mv -n "$src" "$dst"
  if [ -e "$src" ]
  then
    echo "### ERROR: mv failed for $src -- aborting ###" | tee -a "$LOG"
    exit 1
  fi
  echo "moved $src -> $dst" | tee -a "$LOG"
done

# ---- PHASE 1: build L1-L4 binaries ----
echo "---- PHASE 1: build ----" | tee -a "$LOG"
for L in 1 2 3 4
do
  BIN=glue_therm_flow_L${L}_cu129_claude.o
  if [ -f "$BIN" ] && [ "$BIN" -nt "$SRC" ]
  then
    echo "### L${L}: $BIN up-to-date, skip ###" | tee -a "$LOG"
    continue
  fi
  echo "### compile L${L} -> $BIN [$(date '+%F %T')] ###" | tee -a "$LOG"
  $NVCC $NVCCBASE $INCLUDES $LDFLAGS -DN_REFINE_CLI=$L -o "$BIN" "$SRC" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]
  then
    echo "### ERROR: build failed for L=$L (rc=$rc) -- aborting ###" | tee -a "$LOG"
    exit "$rc"
  fi
done

# ---- PHASE 2: 30 ensembles, largest number of configs first ----
echo "---- PHASE 2: ensembles ----" | tee -a "$LOG"
ENS_SORTED=$(
  for d in $PAT_STALE $PAT_AT01
  do
    [ -d "$d" ] || continue
    n=$(find "$d" -maxdepth 1 -name 'ckpoint_lat.*' | wc -l)
    echo "$n $d"
  done | sort -rn | awk '{print $2}'
)
NENS=$(printf '%s\n' "$ENS_SORTED" | grep -c . || true)
echo "ensembles: $NENS (expect 30)" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | tee -a "$LOG"
if [ "$NENS" -ne 30 ]
then
  echo "### ERROR: expected 30 ensembles, found $NENS -- aborting ###" | tee -a "$LOG"
  exit 1
fi

run_one () {
  local ens="$1"
  local L nf gsq at tmax tag elog rc
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  at=$(printf '%s' "$ens" | grep -oE 'at[0-9.]+nu0' | sed 's/at//;s/nu0//')
  tmax=$(awk -v g="$gsq" -v l="$L" 'BEGIN{printf "%.6f", 20.0*g/l}')
  tag="Nf${nf}_gsq${gsq}_at${at}_L${L}"
  elog="therm_flow_${tag}_topup_claude.log"
  echo "[start $(date '+%F %T')] $tag  tmax=$tmax" >> "$LOG"
  ./glue_therm_flow_L${L}_cu129_claude.o "$gsq" "$nf" 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/ "$tmax" >> "$elog" 2>&1
  rc=$?
  echo "[done  $(date '+%F %T')] $tag (status $rc, log: $elog)" >> "$LOG"
}
export -f run_one
export KMAX KMIN STRIDE LOG

echo "---- PHASE 2: launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "$ENS_SORTED" | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}

# ---- PHASE 3: coverage check (flowed configs vs configs on disk) ----
echo "---- PHASE 3: coverage (ens | cfg on disk | cfg in h5 | max k in h5 | n_t | tmax) ----" | tee -a "$LOG"
python3 - <<'PY' 2>&1 | tee -a "$LOG"
import glob
import os
import h5py
pats = ['Nf*_gsq*at0.200000nu01.000000mRe0.000000mIm0.000000nt128L[34]_hb0.400000-1.000000',
        'Nf*_gsq*at0.100000nu01.000000mRe0.000000mIm0.000000nt128L[12]_hb1.000000']
ens_list = []
for p in pats:
    ens_list += sorted(glob.glob(p))
nbad = 0
for ens in ens_list:
    ncfg = len([f for f in os.listdir(ens) if f.startswith('ckpoint_lat.')])
    h = 'data_' + ens + '/therm_flow_claude.h5'
    if not os.path.exists(h):
        print(ens, '| NO H5')
        nbad += 1
        continue
    with h5py.File(h, 'r') as F:
        kl = F['klist'][()]
        nt = len(F['tlist'][()])
        tmax = float(F['tmax'][()][0])
    flag = '' if len(kl) == ncfg else '  <-- INCOMPLETE'
    if flag:
        nbad += 1
    print(ens, '|', ncfg, '|', len(kl), '|', int(kl.max()), '|', nt, '|', tmax, flag)
print('ensembles checked:', len(ens_list), ' INCOMPLETE:', nbad)
PY

echo "================ flow top-up DONE $(date) ================" | tee -a "$LOG"
