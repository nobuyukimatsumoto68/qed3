#!/usr/bin/env bash
# tmp_eig_dw_irl_free_L1_claude.sh
# _claude: VALIDATE the IRL+Rayleigh-Ritz low-D_W finder against the DENSE reference, in the FREE field at L1
# (where we already have the exact dense D_W spectrum, eig_dw_spectrum_L1_free_nt128_claude.dat).  Runs
# eig_lanczos_claude.cu with DW_LOWMODES=1, cold gauge; overlays the IRL Ritz values on the dense low modes.
# If they coincide, the method is trusted for the (dense-infeasible) L2/L4 interacting runs.  GPU 1.
# Reads back: eig_dw_irl_free_L1_claude.log + eig_dw_irl_overlay_free_L1_claude.png .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128
unset CONFIG_LAT          # free field (cold gauge)

DEG=8
BETA=2.0
NK=32
NM=96

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dw_irl_free_L1_claude.log
DENSE="eig_dw_spectrum_L1_free_nt128_claude.dat"   # dense free reference (already computed)

echo "### D_W low-mode IRL VALIDATION vs dense  FREE L1 Nt=${NT}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN="eig_arnoldi_L${LREF}_nt${NT}.o"
echo "### BUILD (Arnoldi driver) -> $BIN ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_arnoldi_claude.cu -o "$BIN" 2>>"$LOG" \
  || { echo "### BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

echo "### RUN DW_ARNOLDI free  m(=Nm)=${NM}  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
DW_ARNOLDI=1 CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" 0.0 0.0 -1 ${BETA} ${DEG} ${NK} ${NM} 2>&1 | tee -a "$LOG" \
  || { echo "### RUN FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### OVERLAY plot (Arnoldi vs dense)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
if [ -f "$DENSE" ]; then
  python3 eig_dw_irl_overlay_claude.py "$DENSE" \
    "eig_dw_arnoldi_L${LREF}_nt${NT}_claude.dat" \
    "eig_dw_arnoldi_overlay_free_L1_claude.png" 1.0 2>&1 | tee -a "$LOG"
else
  echo "### dense reference $DENSE MISSING -- run tmp_eig_dw_dov_L1_free_claude.sh first for the overlay ###" | tee -a "$LOG"
  python3 eig_spectrum_scatter_claude.py "eig_dw_arnoldi_L${LREF}_nt${NT}_claude.dat" \
    "eig_dw_arnoldi_free_L1_claude.png" "Arnoldi low D_W (free L1)" 1.0 2>&1 | tee -a "$LOG"
fi

echo "### free L1 validation done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
