#!/bin/bash
# _claude: L4 RETUNE (2026-07-16, NM) -- massless, MIXED force. New Hasenbusch params in
# tuning/includes/hasenbusch_ladder_claude.h: ladder {0.0,0.2,0.4,0.6,1.0} steps {2,2,2,2,4} tau 1.0 gauge
# MG40 (was {0,0.4,1.0}{4,4,4} tau1.2 MG20; gsq6 was stiff acc0.60). 5-frame + gauge x40 -> FRESH dirs.
# DEFAULT = 4-TRAJ CHECK (KMAX=5): seeds the new dir from the old thermalized config (hot start) and runs 4
# new-ladder trajectories per gsq to read dH/acceptance quickly. KMAX=40 bash ... = the full retune (resumes).
#   L4 gsq{2.0,4.0,6.0}, Nf2, massless, -DMIXED_FORCE (~1.35x).
# Both GPUs, gsq6 PRIORITIZED: gpu0={gsq6.0 SOLO, full GPU}; gpu1={gsq2.0 + gsq4.0 MPS-packed}.
# Read dH/acceptance from the "# dH : .. is_accept : .. rate : .." lines in the tee'd logs.
#
# Run detached:  nohup bash run_retune_L4_claude.sh > run_retune_L4_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/tuning
cd "$SRCDIR" || exit 1
source ../../env.sh

NF=2
NU0=1.0
KMAX=${KMAX:-5}   # DEFAULT = 4-traj check (seed at cfg1 -> k=2..4). KMAX=40 bash ... = full retune (resumes).
KRNG=1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_claude.cu

# ---- MPS daemon ---------------------------------------------------------------------------------
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "MPS daemon: already running"
else
  echo "MPS daemon not up -- starting nvidia-cuda-mps-control -d"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "ERROR: MPS daemon failed to start -- aborting"; exit 1; }

# ---- build L4 (mixed force) ---------------------------------------------------------------------
OUT=hmc_retune_L4.out
echo "===== build $OUT (LREF=4 KMAX=$KMAX -DMIXED_FORCE)  [$(date +%F_%H:%M:%S)] ====="
$NVCC $NVCCFLAGS $INCLUDES -DLREF=4 -DKMAX=$KMAX -DKRNG=$KRNG -DMIXED_FORCE "$SRC" $LDFLAGS -o "$OUT" 2>&1 | tee build_retune_L4_claude.log
test -f "$OUT" || { echo "BUILD FAILED"; exit 1; }

# ---- launch (massless: mass_re=mass_im=0) -------------------------------------------------------
run_one () {   # gpu gsq
  local gpu=$1 gsq=$2
  echo "### [gpu$gpu L4 gsq$gsq] start [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$gpu ./$OUT $gsq $NF $NU0 0.0 0.0 \
    2>&1 | tee run_retune_L4_gsq${gsq}_claude.log
  echo "### [gpu$gpu L4 gsq$gsq] done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
}
# ---- seed the NEW-ladder dir from the OLD run's latest config (hot start; NO re-thermalization needed --
# the Hasenbusch ladder is an exact rewriting of the same determinant, so the old equilibrium config is valid).
# Non-destructive: only copies into a FRESH new dir (skips if it already has configs); no rm/overwrite.
OLDHB="hb0.400000-1.000000"
NEWHB="hb0.200000-0.400000-0.600000-1.000000"
seed () {   # gsq
  local gsq=$1
  local gf=$(printf "%.6f" "$gsq")
  local base="Nf2_gsq${gf}at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4"
  local old="${base}_${OLDHB}"
  local new="${base}_${NEWHB}"
  if ls "$new"/ckpoint_lat.* >/dev/null 2>&1; then
    echo "seed gsq$gsq: $new already has configs -- skip (driver resumes)"
    return
  fi
  local idx=""
  for k in $(ls "$old"/ckpoint_lat.* 2>/dev/null | sed 's/.*\.//' | sort -rn); do
    [ -f "$old/ckpoint_rng.$k" ] && { idx=$k; break; }
  done
  if [ -z "$idx" ]; then echo "seed gsq$gsq: no old lat+rng pair -- COLD start"; return; fi
  # copy the old latest config into the new dir at a FIXED index 1 (fresh numbering; driver resumes from 1,
  # so k=2..KMAX = KMAX-1 new-ladder trajectories -- KMAX=5 -> 4 traj uniform across gsq).
  mkdir -p "$new"
  cp "$old/ckpoint_lat.$idx" "$new/ckpoint_lat.1"
  cp "$old/ckpoint_rng.$idx" "$new/ckpoint_rng.1"
  echo "seed gsq$gsq: hot start from old cfg $idx -> $new/ckpoint_{lat,rng}.1"
}
for g in 2.0 4.0 6.0; do seed "$g"; done

echo "===== launch: gpu0={gsq6.0 SOLO, priority}, gpu1={gsq2.0 + gsq4.0 MPS}  [$(date +%F_%H:%M:%S)] ====="
run_one 0 6.0 &   # gpu0: gsq6 SOLO (stiff coupling, prioritized -- full GPU)
run_one 1 2.0 &   # gpu1 slot A (MPS-packed with gsq4)
run_one 1 4.0 &   # gpu1 slot B
wait
echo "===== L4 retune finished  [$(date +%F_%H:%M:%S)] ====="
