#!/usr/bin/env bash
# run_conn_s2_L2_local_claude.sh -- LOCAL (barracuda22) fill of the L2 at=0.2 stride-2 connected Y_lm tower
# REMAINDER (FNAL-inherited; FNAL allocation consumed 2026-09-01).  Requested by qed3-b7 (2026-10-07).
# Derived from run_conn_s2_local_claude.sh (L3 / L1 at0.1 runner, left untouched).  jj block-avg: arXiv:0804.1501.
#
# Remainder at 2026-10-07 (h0, 400 per residue class k mod 10, total 1717):
#   Nf2 g1.0 cls5 400 | Nf2 g3.0 cls7 68 | Nf4 g1.0 cls7 400 | Nf4 g3.0 cls9 320 | Nf6 g3.0 cls3 328, cls5 153, cls7 48
#
# Differences vs the L3/L1 runner:
#   * TOOLCHAIN (2026-09-30 upgrade): CUDA 12.6 is gone -> /usr/local/cuda-12.9/bin/nvcc -ccbin g++-12, and
#     cuda-12.9 lib64 first on LD_LIBRARY_PATH at run time.
#   * Work list = only the MISSING k's, found at launch by scanning the h5 dir, split into chunks of CHUNK
#     configs (driver --kmin first --kmax last+1 --stride 10).  Only 7 (ensemble, class) units are unfinished
#     and 4 of them are 300-400 configs, so whole units would load-balance badly over 4 workers; chunks fix that.
#     The driver's cheap complete-gate pre-skip makes any already-done k inside a chunk cost ~nothing.
#   * Both GPUs by default (GPU_LIST="0 1") x MPS NPACK=2 = 4 workers.
# Knobs: TPONLY=1 (-DTP_ONLY: s3 + scalar only, no s1/s2, ~2x faster; h5 key A_MIN=3), GPU_LIST, CHUNK (default 25),
#        FORCE_BUILD=1.
# RESUMABLE (complete-gated). NO rm. Stop = kill -- -<PGID of this script>.
#
# Run detached:
#   TPONLY=1 nohup bash run_conn_s2_L2_local_claude.sh > run_conn_s2_L2_local_claude.log 2>&1 &
# Per-worker logs: conn_s2L2_w<wid>_claude.log (this dir).
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:${LD_LIBRARY_PATH:-}

read -ra GPUS <<< "${GPU_LIST:-0 1}"
NGPU=${#GPUS[@]}
NPACK=2
NWORKERS=$(( NGPU * NPACK ))
STRIDE=10
NHITS=1
CHUNK="${CHUNK:-25}"
L=2
AT=0.200000
GLOB='Nf*at0.200000*nt128L2_hb*'
if [ -n "${TPONLY:-}" ]
then
  BINSUF="_tponly"
  DEFEXTRA="-DTP_ONLY"
else
  BINSUF=""
  DEFEXTRA=""
fi

NVCC=/usr/local/cuda-12.9/bin/nvcc
NVCCBASE="-ccbin g++-12 -arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SRC_CONN=jj_local_ylm_scalar_conn_stoch_claude.cu
BIN="jj_local_ylm_scalar_conn_stoch_L${L}${BINSUF}.o"

# ---- build the L2 binary (always rebuilt if missing or FORCE_BUILD=1; else only if the driver source is newer) ----
if [ -n "${FORCE_BUILD:-}" ] || [ ! -f "$BIN" ] || [ "$SRC_CONN" -nt "$BIN" ]
then
  echo "### compile conn L${L} (-DN_REFINE_CLI=${L} $DEFEXTRA) -> $BIN  [$(date +%F_%H:%M:%S)] ###"
  $NVCC $NVCCBASE -DN_REFINE_CLI=${L} $DEFEXTRA $INCLUDES $LDFLAGS "$SRC_CONN" -o "$BIN" \
    || { echo "### conn L${L} BUILD FAILED ###"; exit 1; }
else
  echo "### $BIN up-to-date, skip (FORCE_BUILD=1 to rebuild) ###"
fi
echo "### build OK ###"

# ---- ensure MPS daemon up ----
if pgrep -f nvidia-cuda-mps-control >/dev/null
then
  echo "### MPS daemon: already running ###"
else
  echo "### MPS daemon not up -- starting nvidia-cuda-mps-control -d ###"
  nvidia-cuda-mps-control -d
  for i in 1 2 3 4 5
  do
    pgrep -f nvidia-cuda-mps-control >/dev/null && break
    sleep 1
  done
fi
pgrep -f nvidia-cuda-mps-control >/dev/null \
  || { echo "### ERROR: MPS daemon failed to start -- aborting ###"; exit 1; }

# ---- build the chunked work list: "ens|kmin|kmax" for every run of CHUNK missing k's ----
UNITS=()
# Priority order from the organizer qed3-b7 (2026-10-07): Nf2 g1.0 > Nf4 g1.0 > Nf6 g3.0 > Nf4 g3.0 > Nf2 g3.0,
# then anything else the glob finds (complete ensembles contribute no chunks).
PRIO='Nf2_gsq1.000000 Nf4_gsq1.000000 Nf6_gsq3.000000 Nf4_gsq3.000000 Nf2_gsq3.000000'
ENS_ORDER=()
for pfx in $PRIO
do
  for d in ${pfx}at0.200000*nt128L2_hb*
  do
    [ -d "$d" ] && ENS_ORDER+=("$d")
  done
done
for d in $GLOB
do
  case " ${ENS_ORDER[*]} " in *" $d "*) continue;; esac
  ENS_ORDER+=("$d")
done
for ens in "${ENS_ORDER[@]}"
do
  [ -d "$ens" ] || continue
  case "$ens" in *_vmRe*) continue;; esac
  case "$ens" in *mRe0.000000*) ;; *) continue;; esac
  last=$(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n | tail -1)
  [ -n "$last" ] || continue
  dd="data_${ens}_vmRe0.000000vmIm0.000000/corr_ylm_conn_t00_nhits1_s1"
  for off in 1 3 5 7 9
  do
    missing=$(
      for (( k=off; k<=last; k+=STRIDE ))
      do
        [ -f "$dd/corr.$k.h0.h5" ] || echo "$k"
      done
    )
    [ -n "$missing" ] || continue
    chunks=$(printf '%s\n' "$missing" | awk -v c="$CHUNK" '
      { if (n == 0) first = $1
        lastk = $1
        n++
        if (n == c) { print first, lastk + 1
                      n = 0 } }
      END { if (n > 0) print first, lastk + 1 }')
    while read -r kmin kmax
    do
      UNITS+=("$ens|$kmin|$kmax")
    done <<< "$chunks"
  done
done
echo "### START conn-s2 L2 at=$AT: GPUs=(${GPUS[*]}) x MPS ${NPACK}-pack = ${NWORKERS} workers, ${#UNITS[@]} chunks of <=${CHUNK} cfg, BIN=$BIN  [$(date +%F_%H:%M:%S)] ###"
for u in "${UNITS[@]}"
do
  echo "  unit $u"
done

# ---- workers: round-robin over chunks ----
for (( w=0; w<NWORKERS; w++ ))
do
  (
    gpu=${GPUS[ w % NGPU ]}
    LOG="conn_s2L${L}_w${w}_claude.log"
    idx=$w
    while [ "$idx" -lt "${#UNITS[@]}" ]
    do
      u="${UNITS[$idx]}"
      ens="${u%%|*}"
      rest="${u#*|}"
      kmin="${rest%%|*}"
      kmax="${rest#*|}"
      nf=$(printf '%s' "$ens" | grep -oE '^Nf[0-9]+' | sed 's/Nf//')
      gsq=$(printf '%s' "$ens" | grep -oE 'gsq[0-9.]+at' | sed 's/gsq//;s/at//')
      {
        echo "### CONN-S2 $ens L${L} kmin=$kmin kmax=$kmax stride=$STRIDE at=$AT gpu=$gpu  [$(date +%F_%H:%M:%S)] ###"
        CUDA_VISIBLE_DEVICES=$gpu ./"$BIN" --gsq "$gsq" --Nf "$nf" --nu0 1.0 --at "$AT" \
          --ens-dir "$ens/" --kmin "$kmin" --stride "$STRIDE" --kmax "$kmax" \
          --nhits "$NHITS" --t0 0 --spin-dilution
        echo "### CONN-S2 $ens kmin=$kmin done (status $?)  [$(date +%F_%H:%M:%S)] ###"
      } >> "$LOG" 2>&1
      idx=$(( idx + NWORKERS ))
    done
  ) &
done
wait
echo "### ALL conn-s2 L2 done  [$(date +%F_%H:%M:%S)] ###"
