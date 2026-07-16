#!/usr/bin/env bash
# tmp_eig_validate_claude.sh
# _claude: VALIDATE the pre-existing IRL (includes/lanczos_claude.h) against the dense geev reference, on the
# free field at L1 Nt32 (N=768).  Builds eig_lanczos_claude.cu and eig_wmass_val_claude.cu, runs both massless
# (mass=0), then diffs the lowest Nk eigenvalues of A=(D_ov+m)^dag(D_ov+m).  Single process on GPU 0.
# Reads back: eig_validate_claude.log  and the two .dat files + the PASS/FAIL line.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=32
NK=1          # free field is 8-fold DEGENERATE at the lowest level (dense: eight modes at 0.96084); single-
NM=16         # vector Lanczos resolves ONE copy cleanly, so Nk=1 is the bulletproof machinery check.  Multi-
DEG=12        # mode needs a NON-degenerate spectrum (interacting config).  EVEN degree (lanczos-2.pdf).
ALPHA=-1      # auto = 2+|m| (overlap singular max ~2)
BETA=1.15     # above the lowest singular level (0.98), below the next (1.285)

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_validate_claude.log
# (log appends; '### ' headers + timestamps delimit each run.  Rename the old log if you want a clean one.)

echo "### VALIDATE IRL vs dense geev   L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

# ---- build IRL driver ----
BIN_IRL="eig_lanczos_L${LREF}.o"
echo "### BUILD IRL -> $BIN_IRL ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_lanczos_claude.cu -o "$BIN_IRL" 2>>"$LOG" \
  || { echo "### IRL BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

# ---- build dense reference ----
BIN_DEN="eig_wmass_val_L${LREF}.o"
echo "### BUILD dense -> $BIN_DEN ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DEN" 2>>"$LOG" \
  || { echo "### dense BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

# ---- run dense reference (massless) ----
echo "### RUN dense (mass=0)  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DEN" 0.0 0.0 2>&1 | tee -a "$LOG" \
  || { echo "### dense RUN FAILED ###" | tee -a "$LOG"; exit 1; }

# ---- run IRL (massless), auto Chebyshev window ----
echo "### RUN IRL (mass=0, alpha=${ALPHA} beta=${BETA} Nk=${NK} Nm=${NM} deg=${DEG})  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_IRL" 0.0 0.0 ${ALPHA} ${BETA} ${DEG} ${NK} ${NM} 2>&1 | tee -a "$LOG" \
  || { echo "### IRL RUN FAILED ###" | tee -a "$LOG"; exit 1; }

# ---- diff the lowest Nk ----
echo "### COMPARE lowest ${NK}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 compare_eig_validate_claude.py \
  "eig_lanczos_L${LREF}_nt${NT}_claude.dat" \
  "eig_wmass_L${LREF}_nt${NT}_claude.dat" \
  ${NK} 2>&1 | tee -a "$LOG"

echo "### validate done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
