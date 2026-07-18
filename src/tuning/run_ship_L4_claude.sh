#!/bin/bash
# _claude: SHIPPED-L4 4-traj SANITY (2026-07-16). Confirms the shipped L4 config -- {0,0.4,1.0} steps {4,4,4}
# tau 0.8 gauge MG100 (from tuning/includes/hasenbusch_ladder/steps/tau/mg) -- gives healthy dH/acceptance at
# the stiff gsq6.0 before full production. Mixed force. Just builds the tuning driver with the current header.
#
# ISOLATED + hot-start: writes a FRESH dir test_ship_L4_gsq<g>/ via the OUTDIR env, SEEDED from that gsq's
# latest existing _hb0.4-1.0 config copied to index 1 -> driver resumes from 1, runs k=2..5 = 4 trajectories.
# (Gauge config is ladder/MG/tau-independent -> valid hot start, no re-thermalization.)
#
# HANDOFF -- the USER runs it:  nohup bash run_ship_L4_claude.sh > run_ship_L4_claude.log 2>&1 &
#   (GSQ env picks the coupling, default 6.0 = stiffest; GPU env, default 1. GSQ=2.0 bash ... for another.)
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/tuning
cd "$SRCDIR" || exit 1
source ../../env.sh

NF=2
NU0=1.0
GSQ=${GSQ:-6.0}
GPU=${GPU:-1}
KMAX=5            # traj = KMAX-1 = 4 (seed at index 1)
KRNG=1
SRCHB="hb0.400000-1.000000"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_claude.cu

# ---- MPS daemon (single proc, but keep consistent) ----------------------------------------------
if pgrep -f nvidia-cuda-mps-control >/dev/null
then echo "MPS daemon: already running"
else
  echo "MPS daemon not up -- starting nvidia-cuda-mps-control -d"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
fi

# ---- build (shipped header, mixed force) --------------------------------------------------------
BIN=hmc_ship_L4.out
echo "===== build $BIN (LREF=4 KMAX=$KMAX -DMIXED_FORCE)  [$(date +%F_%H:%M:%S)] ====="
$NVCC $NVCCFLAGS $INCLUDES -DLREF=4 -DKMAX=$KMAX -DKRNG=$KRNG -DMIXED_FORCE "$SRC" $LDFLAGS -o "$BIN" 2>&1 | tee build_ship_L4_claude.log
test -f "$BIN" || { echo "BUILD FAILED"; exit 1; }

# ---- seed a fresh OUTDIR from gsq's latest config, then run 4 traj -------------------------------
GF=$(printf "%.6f" "$GSQ")
SRCDIR_CFG="Nf2_gsq${GF}at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_${SRCHB}"
OUT="test_ship_L4_gsq${GF}/"
if ! ls "$OUT"ckpoint_lat.* >/dev/null 2>&1; then
  IDX=""
  for k in $(ls "$SRCDIR_CFG"/ckpoint_lat.* 2>/dev/null | sed 's/.*\.//' | sort -rn); do
    [ -f "$SRCDIR_CFG/ckpoint_rng.$k" ] && { IDX=$k; break; }
  done
  [ -z "$IDX" ] && { echo "no src lat+rng in $SRCDIR_CFG -- abort"; exit 1; }
  mkdir -p "$OUT"
  cp "$SRCDIR_CFG/ckpoint_lat.$IDX" "${OUT}ckpoint_lat.1"
  cp "$SRCDIR_CFG/ckpoint_rng.$IDX" "${OUT}ckpoint_rng.1"
  echo "seed: hot start from $SRCDIR_CFG cfg $IDX -> ${OUT}ckpoint_{lat,rng}.1"
else
  echo "seed: $OUT already seeded -- resume"
fi

echo "===== run L4 gsq$GSQ SHIPPED config (tau0.8 MG100), 4 traj on gpu$GPU  [$(date +%F_%H:%M:%S)] ====="
CUDA_VISIBLE_DEVICES=$GPU OUTDIR="$OUT" ./$BIN $GSQ $NF $NU0 0.0 0.0 \
  2>&1 | tee run_ship_L4_gsq${GSQ}_claude.log
echo "===== done  [$(date +%F_%H:%M:%S)] ====="
