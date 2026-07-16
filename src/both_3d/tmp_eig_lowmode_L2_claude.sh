#!/usr/bin/env bash
# tmp_eig_lowmode_L2_claude.sh
# _claude: IRL low-mode run at L2 (N=10752, dense geev already INFEASIBLE -> self-validating IRL only).
# Intermediate step toward the L4 deflation/LMA target.  Runs on the interacting config with the SMALLEST
# Wilson lambda_min (ckpoint_lat.5, lambda_min=0.0485 = strongest near-zero tail) so the low modes deflation
# would exploit are most pronounced.  The driver prints per-mode relative residual + PASS/FAIL and writes the
# low spectrum to eig_lanczos_L2_nt128_claude.dat.  Single process on GPU 1.  Reads back:
# eig_lowmode_L2_claude.log .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=2
NT=128
CONFIG="Nf2_gsq8.000000at0.200000nu01.000000nt128L2/ckpoint_lat.740"  # THERMALIZED (traj 740; ens runs to
              # 1601).  Its Wilson lambda_min=0.0643 sits ABOVE the frozen window edge 0.06, so the near-zero
              # mode is COVERED by the Zolotarev fit (the earlier k=5 was traj-5 NON-thermalized, lambda_min
              # 0.0485 < 0.06 -> below window -> near-zero overlap mode distorted -> IRL found only bulk ~3.5).

NK=4          # lowest 4 near-zero modes (tight beta=0.5 band -> few modes; k small for clean convergence)
NM=16         # Krylov dim (4x k)
DEG=8         # EVEN.  deg=12 gave ~2.7e5 gain at mu=0.05 (excessive dynamic range -> Lanczos stress); deg=8
              # gives gentler gain (~20-30x at beta=0.5 vs bulk <=1).  See cheby_filter_claude.png.
ALPHA=-1      # auto = 2+|m| (overlap singular max ~2; header squares -> alpha^2=4 = max eigenvalue of Dov^2)
BETA=0.5      # wanted band = mu < beta^2 = 0.25 (|z| < 0.5): focus on the near-zero TAIL, exclude bulk-low.
              # NB gentler gain here (deg=8 -> ~20-30x vs bulk<=1; see cheby_filter_claude.png) since q(0)
              # only just below -1; still enough separation with full re-orth.

export CONFIG_LAT="$CONFIG"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_lowmode_L2_claude.log
# (log appends; '### ' headers + timestamps delimit each run.)

echo "### IRL low-mode L2 (self-validating)   Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  config=${CONFIG}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN_IRL="eig_lanczos_L${LREF}_nt${NT}.o"
echo "### BUILD IRL -> $BIN_IRL ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_lanczos_claude.cu -o "$BIN_IRL" 2>>"$LOG" \
  || { echo "### IRL BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

echo "### RUN IRL (mass=0, alpha=${ALPHA} beta=${BETA} Nk=${NK} Nm=${NM} deg=${DEG})  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_IRL" 0.0 0.0 ${ALPHA} ${BETA} ${DEG} ${NK} ${NM} 2>&1 | tee -a "$LOG" \
  || { echo "### IRL RUN FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### L2 low spectrum + residuals (eig_lanczos_L${LREF}_nt${NT}_claude.dat) ###" | tee -a "$LOG"
grep -v '^#' "eig_lanczos_L${LREF}_nt${NT}_claude.dat" | tee -a "$LOG"

echo "### L2 low-mode done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
