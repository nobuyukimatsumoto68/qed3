#!/usr/bin/env bash
# run_cond_L4_match_claude.sh -- measure the scalar condensate on the L4 massive ensembles at the
# SAME physical masses as L1/L2 (m = 0.1, 0.2, 0.3, 0.4) so cond_vs_mass has matched-mass points at
# every L. Reuses jj_scalar_condensate_eo_L4.o (e/o+spin dilution).  m=0.1 already has 3 cfg
# (complete-gated -> resumes/tops up); 0.2/0.3/0.4 fresh.  stride 4 (L4 streams ~59 cfg -> ~15 pts),
# kmin 20 (therm cut), nhits 1.  NO MPS -- one process at a time on GPU 0.  No rm anywhere.
#
# Run:  nohup bash run_cond_L4_match_claude.sh > run_cond_L4_match_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
GPU="${GPU:-0}"
STRIDE=4
KMIN=20
NHITS=1
BIN=jj_scalar_condensate_eo_L4.o

# ensemble dir (col-1 basename) per matched mass -- hb tag is the resc-shifted one for L4 gsq6
declare -A ENS=(
  [0.100000]="Nf2_gsq6.000000at0.200000nu01.000000mRe0.100000mIm0.000000nt128L4_hb0.425902-1.025902"
  [0.200000]="Nf2_gsq6.000000at0.200000nu01.000000mRe0.200000mIm0.000000nt128L4_hb0.451804-1.051804"
  [0.300000]="Nf2_gsq6.000000at0.200000nu01.000000mRe0.300000mIm0.000000nt128L4_hb0.477706-1.077706"
  [0.400000]="Nf2_gsq6.000000at0.200000nu01.000000mRe0.400000mIm0.000000nt128L4_hb0.503608-1.103608"
)

for m in 0.100000 0.200000 0.300000 0.400000
do
  ens="${ENS[$m]}"
  LOG="cond_L4_match_mRe${m}_claude.log"
  if [ ! -d "$ens" ]
  then
    echo "### SKIP m=$m: dir absent ($ens)  [$(date +%T)] ###" | tee -a "$LOG"
    continue
  fi
  ks=$(ls "$ens"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
  last=$(echo "$ks" | tail -1)
  {
    echo "### COND L4 m=$m  ens=$ens  kmin=$KMIN stride=$STRIDE kmax=$((last+1))  [$(date +%T)] ###"
    CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" --ens-dir "$ens/" --kmin "$KMIN" --stride "$STRIDE" \
      --kmax "$((last+1))" --nhits "$NHITS" --mass-re "$m"
    echo "### done m=$m (status $?)  [$(date +%T)] ###"
  } >> "$LOG" 2>&1
  echo "done m=$m"
done
echo "### ALL L4 matched-mass condensate done  [$(date +%T)] ###"
