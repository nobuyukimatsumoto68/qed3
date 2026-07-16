#!/bin/bash
# _claude: 4-traj ACCEPTANCE TEST, L=4 gsq6.0 (largest), two 5-frame ladder CONFIGS head-to-head, MPS-packed
# (2026-07-16, NM). Both gauge MG40 (hasenbusch_mg L4->40), tau 1.0, mixed force. Compile-time -DLADDER_B picks:
#   A (default): {0,0.2,0.4,0.6,1.0}  steps {2,2,2,2,4}  (heaviest x2)
#   B (-DLADDER_B): {0,0.2,0.4,0.55,1.0} steps {2,2,2,2,6}  (heaviest x3, 4th rung 0.55)
#
# ISOLATED + hot-start: each config writes a FRESH dir (test_5f_A_L4_gsq6/ , test_5f_B_L4_gsq6/) SEEDED from
# gsq6's latest existing _hb0.4-1.0 config copied to index 1 -> driver resumes from 1, runs k=2..5 = 4 traj.
# (Gauge config is ladder/MG-independent -> valid hot start, no re-thermalization.)
#
# MPS: A and B packed together on ONE GPU (GPU env, default 1). Read dH/acceptance from the tee'd logs.
#
# THIS IS A HANDOFF -- the USER runs it:  nohup bash run_test_mg40_claude.sh > run_test_mg40_claude.log 2>&1 &
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/tuning
cd "$SRCDIR" || exit 1
source ../../env.sh

NF=2
NU0=1.0
GSQ=6.0
KMAX=5            # traj = KMAX-1 = 4 (seed at index 1)
KRNG=1
GPU=${GPU:-1}     # both configs MPS-packed on this GPU
SRCHB="hb0.400000-1.000000"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=hmc_hasenbusch_claude.cu

# ---- MPS daemon ---------------------------------------------------------------------------------
if pgrep -f nvidia-cuda-mps-control >/dev/null
then echo "MPS daemon: already running"
else
  echo "MPS daemon not up -- starting nvidia-cuda-mps-control -d"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5; do pgrep -f nvidia-cuda-mps-control >/dev/null && break; sleep 1; done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null || { echo "ERROR: MPS daemon failed to start"; exit 1; }

# ---- build the two config binaries (A default, B = -DLADDER_B); both LREF=4 mixed force -----------
build () {   # cfg extra
  local cfg=$1 extra=$2
  local out=hmc_5f_${cfg}_L4.out
  echo "===== build $out (LREF=4 KMAX=$KMAX -DMIXED_FORCE $extra)  [$(date +%F_%H:%M:%S)] ====="
  $NVCC $NVCCFLAGS $INCLUDES -DLREF=4 -DKMAX=$KMAX -DKRNG=$KRNG -DMIXED_FORCE $extra "$SRC" $LDFLAGS -o "$out" 2>&1 | tee build_5f_${cfg}_claude.log
  test -f "$out" || { echo "BUILD FAILED: $out"; exit 1; }
}
build A ""
build B "-DLADDER_B"

# ---- seed a fresh OUTDIR test dir from gsq6's latest config, then run 4 traj ----------------------
seed_and_run () {   # cfg
  local cfg=$1
  local gf=$(printf "%.6f" "$GSQ")
  local src="Nf2_gsq${gf}at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_${SRCHB}"
  local out="test_5f_${cfg}_L4_gsq${gf}/"
  if ! ls "$out"ckpoint_lat.* >/dev/null 2>&1; then
    local idx=""
    for k in $(ls "$src"/ckpoint_lat.* 2>/dev/null | sed 's/.*\.//' | sort -rn); do
      [ -f "$src/ckpoint_rng.$k" ] && { idx=$k; break; }
    done
    if [ -z "$idx" ]; then echo "seed $cfg: no src lat+rng in $src -- SKIP"; return; fi
    mkdir -p "$out"
    cp "$src/ckpoint_lat.$idx" "${out}ckpoint_lat.1"
    cp "$src/ckpoint_rng.$idx" "${out}ckpoint_rng.1"
    echo "seed $cfg: hot start from $src cfg $idx -> ${out}ckpoint_{lat,rng}.1"
  else
    echo "seed $cfg: $out already seeded -- resume"
  fi
  echo "### [gpu$GPU cfg$cfg gsq$GSQ] start [$(date +%F_%H:%M:%S)] ###"
  CUDA_VISIBLE_DEVICES=$GPU OUTDIR="$out" ./hmc_5f_${cfg}_L4.out $GSQ $NF $NU0 0.0 0.0 \
    2>&1 | tee run_test_5f_${cfg}_L4_gsq${GSQ}_claude.log
  echo "### [gpu$GPU cfg$cfg gsq$GSQ] done (status ${PIPESTATUS[0]}) [$(date +%F_%H:%M:%S)] ###"
}

echo "===== launch: gpu$GPU = {config A + config B MPS-packed} at gsq$GSQ  [$(date +%F_%H:%M:%S)] ====="
seed_and_run A &
seed_and_run B &
wait
echo "===== 5-frame A/B 4-traj test finished  [$(date +%F_%H:%M:%S)] ====="
