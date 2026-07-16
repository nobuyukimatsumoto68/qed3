#!/usr/bin/env bash
# tmp_eig_validate_interacting_claude.sh
# _claude: MULTI-MODE IRL validation on a real (INTERACTING) L1 config, whose low overlap spectrum is
# NON-degenerate (unlike the free field, which is 8-fold degenerate and only supports an Nk=1 check).
# NO DENSE geev here -- the IRL is SELF-VALIDATING: eig_lanczos recomputes, per returned mode, the relative
# residual ||A u - lam u|| / (lam ||u||), which is the DEFINITION of an eigenpair (PASS if all < 1e-7).  That
# replaces the O(N^2) dense reference (the cheap free-field N=768 dense already anchors that the filter finds
# the SMALLEST eigenvalue).  Builds eig_lanczos_claude.cu at L1 Nt128 (N=3072), points it at the config via
# CONFIG_LAT, runs the IRL.  Single process on GPU 1.  Reads back: eig_validate_interacting_claude.log.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128
CONFIG="Nf2_gsq8.000000at0.200000nu01.000000nt128L1/ckpoint_lat.2000"

NK=8          # interacting low modes are non-degenerate -> multi-mode Lanczos resolves them
NM=48
DEG=12        # EVEN (lanczos-2.pdf)
ALPHA=-1      # auto = 2+|m| (overlap singular max ~2)
BETA=1.8      # wanted band = singular < 1.8 (mu < 3.24).  This config is GAPPED: smallest overlap singular
              # ~1.49 (mu~2.227), NO near-zero modes at L1 gsq8.  beta=1.5 (mu<2.25) bracketed only the
              # lowest 4 (2.227-2.244) -> 4/8; the 5th sits just above 2.25.  beta=1.8 brackets >=8 with
              # good filter contrast (lowest 8 at q~-3.6, bulk >3.4 at q>-1 suppressed).

export CONFIG_LAT="$CONFIG"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_validate_interacting_claude.log
# (log appends; '### ' headers + timestamps delimit each run.)

echo "### INTERACTING IRL validation (self-validating, NO dense)   L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  config=${CONFIG}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

# ---- build IRL driver ----
BIN_IRL="eig_lanczos_L${LREF}_nt${NT}.o"
echo "### BUILD IRL -> $BIN_IRL ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_lanczos_claude.cu -o "$BIN_IRL" 2>>"$LOG" \
  || { echo "### IRL BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

# ---- run IRL (massless); the driver prints per-mode residual + PASS/FAIL ----
echo "### RUN IRL (mass=0, alpha=${ALPHA} beta=${BETA} Nk=${NK} Nm=${NM} deg=${DEG})  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_IRL" 0.0 0.0 ${ALPHA} ${BETA} ${DEG} ${NK} ${NM} 2>&1 | tee -a "$LOG" \
  || { echo "### IRL RUN FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### low spectrum + residuals (eig_lanczos_L${LREF}_nt${NT}_claude.dat) ###" | tee -a "$LOG"
grep -v '^#' "eig_lanczos_L${LREF}_nt${NT}_claude.dat" | tee -a "$LOG"

echo "### interacting validate done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
