#!/usr/bin/env bash
# tmp_flow_weighted_claude.sh -- LINK-WEIGHTED spatial Wilson flow (NM paper note, 2026-10-26),
# first pass: the 12 massless Nf2 at0.2 ensembles (L1-L4, 3 couplings each).
# Flow:  d theta_l/d tau_hat = -(ell_l^2/(2 A_l)) d E_hat_t/d theta_l,  E_hat_t = sum_f (abar^2/A_f) theta_f^2/2,
#        abar = base.mean_ell, Luscher RK3 (arXiv:1006.4518 App. C), 200 steps to a per-L
#        tmax = R^2/(2 abar^2) (R=1; flow radius 2R), eps = tmax/200, 201 saved points.
#        tmax: L1 0.408, L2 1.432, L3 3.158, L4 5.575; eps: L1 0.0020 ... L4 0.0279.
# Driver glue_flow_weighted_claude.cu + includes/flow_weighted_claude.h. Output data_<ens>/flow_weighted_claude.h5
# (new file name; the old therm_flow_claude.h5 are untouched). Resumable per config.
# CPU-only, single-threaded workers, one ensemble per worker.
# TOOLCHAIN (2026-09-30): /usr/local/cuda-12.9/bin/nvcc -ccbin g++-12; cuda-12.9 lib64 first.
#
# Run (detached):  nohup bash tmp_flow_weighted_claude.sh > /dev/null 2>&1 &
# Log:             tmp_flow_weighted_claude.log  (per-ensemble: flow_weighted_<tag>_claude.log)
# Workers:         NWORK=12 by default

set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1

export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1

NWORK=${NWORK:-12}
KMAX=100000
KMIN=1
STRIDE=1
# TMAX=4.0   # first draft: uniform
TMAX=-1      # <= 0 -> driver uses the per-L value R^2/(2 abar^2)

NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCBASE="-ccbin g++-12 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lm'
SRC=glue_flow_weighted_claude.cu
HDR=includes/flow_weighted_claude.h

ENS_LIST="
Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000
Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000
Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000
Nf2_gsq4.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000
Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
Nf2_gsq4.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
Nf2_gsq6.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
"

LOG=tmp_flow_weighted_claude.log
echo "================ weighted flow START $(date) ================" | tee -a "$LOG"

# ---- PHASE 1: build L1-L4 ----
echo "---- PHASE 1: build ----" | tee -a "$LOG"
for L in 1 2 3 4
do
  BIN=glue_flow_weighted_L${L}_claude.o
  if [ -f "$BIN" ] && [ "$BIN" -nt "$SRC" ] && [ "$BIN" -nt "$HDR" ]
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

# ---- PHASE 2: smoke test -- first 2 configs of L1 g0.5 and of L4 g2.0 (largest RK3 step).
#      Written into the production h5; resume skips them later. Aborts if non-finite or not decreasing.
echo "---- PHASE 2: smoke (L1 g0.5 and L4 g2.0, k=1..2) ----" | tee -a "$LOG"
SMOKE1=Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000
SMOKE4=Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000
./glue_flow_weighted_L1_claude.o 0.5 2 1.0 2 1 1 "$SMOKE1"/ "$TMAX" 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: smoke run L1 failed (rc=$rc) -- aborting ###" | tee -a "$LOG"
  exit "$rc"
fi
./glue_flow_weighted_L4_claude.o 2.0 2 1.0 2 1 1 "$SMOKE4"/ "$TMAX" 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: smoke run L4 failed (rc=$rc) -- aborting ###" | tee -a "$LOG"
  exit "$rc"
fi
python3 - "$SMOKE1" "$SMOKE4" <<'PY' 2>&1 | tee -a "$LOG"
import sys
import h5py
import numpy as np
ok_all = True
for ens in sys.argv[1:]:
    with h5py.File("data_" + ens + "/flow_weighted_claude.h5", "r") as F:
        t = F["tlist"][()]
        E = F["E"][()]
        k = F["klist"][()]
        abar = float(F["abar"][()][0])
    print("smoke:", ens)
    print("  n_t", len(t), "t_max", t[-1], "eps", t[1] - t[0], "klist", k.tolist(), "abar", abar)
    for i in [0, 1, 10, 50, 100, 200]:
        print("  tau_hat=%.4f  E_hat=%s" % (t[i], np.array2string(E[i], precision=6)))
    finite = bool(np.all(np.isfinite(E)))
    decr = bool(np.all(np.diff(E, axis=0) <= 0))
    print("  finite:", finite, " monotone decreasing in tau_hat:", decr)
    ok_all = ok_all and finite and decr
print("SMOKE_OK" if ok_all else "SMOKE_FAIL")
sys.exit(0 if ok_all else 1)
PY
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: smoke check FAILED (non-finite or non-decreasing E_hat) -- aborting ###" | tee -a "$LOG"
  exit 1
fi

# ---- PHASE 3: the 12 ensembles ----
echo "---- PHASE 3: 12 ensembles, $NWORK workers $(date) ----" | tee -a "$LOG"
run_one () {
  local ens="$1"
  local L gsq tag elog rc
  L=$(printf '%s' "$ens" | grep -oE 'nt128L[0-9]+' | sed 's/nt128L//')
  gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
  tag="Nf2_gsq${gsq}_at0.2_L${L}"
  elog="flow_weighted_${tag}_claude.log"
  echo "[start $(date '+%F %T')] $tag" >> "$LOG"
  ./glue_flow_weighted_L${L}_claude.o "$gsq" 2 1.0 "$KMAX" "$KMIN" "$STRIDE" "$ens"/ "$TMAX" >> "$elog" 2>&1
  rc=$?
  echo "[done  $(date '+%F %T')] $tag (status $rc, log: $elog)" >> "$LOG"
}
export -f run_one
export KMAX KMIN STRIDE TMAX LOG
printf '%s\n' $ENS_LIST | xargs -P "$NWORK" -I{} bash -c 'run_one "$1"' _ {}

# ---- PHASE 4: coverage ----
echo "---- PHASE 4: coverage (ens | cfg on disk | cfg in h5 | max k | n_t | tmax | abar) ----" | tee -a "$LOG"
python3 - $ENS_LIST <<'PY' 2>&1 | tee -a "$LOG"
import os
import sys
import h5py
nbad = 0
for ens in sys.argv[1:]:
    ncfg = len([f for f in os.listdir(ens) if f.startswith("ckpoint_lat.")])
    h = "data_" + ens + "/flow_weighted_claude.h5"
    if not os.path.exists(h):
        print(ens, "| NO H5")
        nbad += 1
        continue
    with h5py.File(h, "r") as F:
        kl = F["klist"][()]
        nt = len(F["tlist"][()])
        tmax = float(F["tmax"][()][0])
        abar = float(F["abar"][()][0])
    flag = "" if len(kl) == ncfg else "  <-- INCOMPLETE"
    if flag:
        nbad += 1
    print(ens, "|", ncfg, "|", len(kl), "|", int(kl.max()), "|", nt, "|", tmax, "|", abar, flag)
print("ensembles checked:", len(sys.argv) - 1, " INCOMPLETE:", nbad)
PY

echo "================ weighted flow DONE $(date) ================" | tee -a "$LOG"
