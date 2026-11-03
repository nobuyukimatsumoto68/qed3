#!/usr/bin/env bash
# run_glue_msm_fill_claude.sh -- FILL the glue_msm_shapes (linear-F shape basis) measurement where it stops before
# the gauge configs (NM 2026-10-12: glue is cheap, extend it to the full k range of the configs/fermion data).
# Missing (checked 2026-10-12 against ckpoint_lat.<k>; every earlier k present, last existing h5 has "complete"):
#   L3 Nf4 g1.5  k 971-999 (29)    L3 Nf4 g3.0  k 965-999 (35)    L3 Nf4 g4.5  k 893-999 (107)
#   L3 Nf6 g1.5  k 956-999 (44)    L3 Nf6 g3.0  k 871-999 (129)   L3 Nf6 g4.5  k 918-999 (82)
#   L4 Nf6 g6.0  k 779-799 (21)                                   total 447 configs
# Uses the EXISTING production binaries glue2_msm_shapes_L{3,4}_claude.o (built 2026-09-14 17:04 from the current
# glue2_msm_shapes_claude.cu: N_FLOW = 1, prefix glue_msm_shapes) -> same code as all existing glue data; NO rebuild.
# Existing outputs are never overwritten: the driver skips every k whose h5 already has the "complete" flag, and each
# job starts at the first missing k.  NO rm, no kill.  CPU-only, single-threaded workers, one ensemble per worker.
# Run from anywhere (cds to src/production; the driver reads ../../geometry/data relative to it):
#   NWORK=4 nohup bash run_glue_msm_fill_claude.sh > run_glue_msm_fill_claude.out 2>&1 &
# All output -> run_glue_msm_fill_claude.log (+ one log per ensemble)
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64:/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib:${LD_LIBRARY_PATH:-}
NWORK="${NWORK:-4}"
LOG=run_glue_msm_fill_claude.log
echo "================ glue_msm_shapes FILL START $(date) NWORK=$NWORK ================" | tee -a "$LOG"

for L in 3 4
do
  if [ ! -x "glue2_msm_shapes_L${L}_claude.o" ]
  then
    echo "ABORT: missing glue2_msm_shapes_L${L}_claude.o" | tee -a "$LOG"
    exit 1
  fi
done

# ensemble | L | gsq | Nf | first missing k | last config k
JOBS=$(cat <<'EOF'
Nf4_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000 3 1.5 4 971 999
Nf4_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000 3 3.0 4 965 999
Nf4_gsq4.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000 3 4.5 4 893 999
Nf6_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000 3 1.5 6 956 999
Nf6_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000 3 3.0 6 871 999
Nf6_gsq4.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L3_hb0.400000-1.000000 3 4.5 6 918 999
Nf6_gsq6.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000 4 6.0 6 779 799
EOF
)

run_one () {
  local ens=$1
  local L=$2
  local gsq=$3
  local nf=$4
  local k0=$5
  local k1=$6
  local elog="glue_msm_fill_Nf${nf}_g${gsq}_L${L}_claude.log"
  local t0
  local t1
  if [ ! -e "$ens/ckpoint_lat.$k1" ]
  then
    echo "ABORT $ens: ckpoint_lat.$k1 missing"
    return
  fi
  t0=$(date +%s)
  # args: gsq Nf nu0 kmax_run kmin stride ens_dir
  ./glue2_msm_shapes_L${L}_claude.o "$gsq" "$nf" 1.0 "$k1" "$k0" 1 "$ens" > "$elog" 2>&1
  t1=$(date +%s)
  local nnew
  nnew=$(grep -c '# read from k' "$elog")
  echo "done Nf${nf} g${gsq} L${L}: k ${k0}-${k1}, measured ${nnew} configs, $((t1 - t0)) s"
}
export -f run_one

printf '%s\n' "$JOBS" | xargs -P "$NWORK" -L 1 bash -c 'run_one "$1" "$2" "$3" "$4" "$5" "$6"' _ 2>&1 | tee -a "$LOG"

echo "---- post-check: glue_msm_shapes h5 vs gauge configs ----" | tee -a "$LOG"
printf '%s\n' "$JOBS" | while read -r ens L gsq nf k0 k1
do
  n=0
  for k in $(seq "$k0" "$k1")
  do
    [ -e "data_$ens/glue_msm_shapes.$k.h5" ] && n=$((n + 1))
  done
  echo "Nf${nf} g${gsq} L${L}: $n of $((k1 - k0 + 1)) filled" | tee -a "$LOG"
done
echo "================ FILL DONE $(date) ================" | tee -a "$LOG"
