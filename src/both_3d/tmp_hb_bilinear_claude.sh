#!/bin/bash
#SBATCH --account=affine.lq2_gpu
#SBATCH --qos=normal
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:20:00
#SBATCH --job-name=hb_bilin_l1
#SBATCH --output=slurm_%x_%j.out
# ---------------------------------------------------------------------------
# _claude HANDOFF (2026-07-11) -- C1 gate for the Hasenbusch external-bra force.
# Builds + runs the L=1 force-vs-FD test (test_hasenbusch_bilinear_l1_claude.cu), which checks
# that OverlapWMass::precalc_grad_bilinear_..._ms + grad reproduce Term B = 2 Re<phi|K|eta>
# (bra != ket) against the central difference of 2 Re<phi|D_ov|eta>. qos=test (<=30 min), affine.
# Run:  sbatch tmp_hb_bilinear_claude.sh     Read back: test_hasenbusch_bilinear_l1_claude.log
# NO rm (nvcc overwrites the .o in place).
# ---------------------------------------------------------------------------
set -u
SRC=/project/qed3/qed3/src/both_3d
CU=test_hasenbusch_bilinear_l1_claude.cu
OBJ=test_hasenbusch_bilinear_l1_claude.o
LOG=$SRC/test_hasenbusch_bilinear_l1_claude.log
ENVSH=/home/nmatsum/env.sh; [ -f "$ENVSH" ] || ENVSH=/lustre2/affine/env.sh

{
echo "######## Hasenbusch bilinear L=1 FD gate  $(date) ########"; hostname; nvidia-smi -L
cd "$SRC" || { echo "no $SRC"; exit 1; }
source "$ENVSH"

echo; echo "== build $OBJ (default reference grad; sm_80) =="
# Extract the exact nvcc line from Makefile.fnal (same idiom as launch_block_build) and run it,
# so we compile ONLY this target with the project flags/includes (no full-tree .d cascade).
LINE=$(make -f Makefile.fnal -n "$OBJ" 2>/dev/null | grep -E "nvcc .*${CU}" | tail -1)
[ -n "$LINE" ] || { echo "ABORT: could not extract the compile line for $OBJ"; exit 1; }
echo "  $LINE" | cut -c1-240
eval "$LINE" || { echo "BUILD FAILED -- ABORT"; exit 1; }
[ -x "$OBJ" ] || { echo "MISSING after build: $OBJ -- ABORT"; exit 1; }
ls -la "$OBJ"

echo; echo "== run =="
OMP_NUM_THREADS=4 ./"$OBJ"
rc=$?
echo; echo "== exit code $rc (0 = PASS) =="
echo "######## DONE ########"
} 2>&1 | tee "$LOG"
