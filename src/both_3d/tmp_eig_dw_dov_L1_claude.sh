#!/usr/bin/env bash
# tmp_eig_dw_dov_L1_claude.sh
# _claude: dense complex spectra of the bare Wilson D_W AND the overlap D_ov for the SAME thermalized L1
# config (gsq8, ckpoint_lat.2105), side by side.  D_ov is a unitary projector of D_W (keeps phase, discards
# |lambda|) -- this shows why the D_W "0.2" mode does NOT give a small |z|.  Two builds of eig_wmass_val_claude.cu
# (-DSPECTRUM_DW, -DSPECTRUM_DOV), each geev at N=3072 (~1 GB, fits), then a 2-panel plot.  GPU 1.
# Reads back: eig_dw_dov_L1_claude.log + eig_dw_dov_compare_claude.png .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128
CONFIG="Nf2_gsq8.000000at0.200000nu01.000000nt128L1/ckpoint_lat.2105"
export CONFIG_LAT="$CONFIG"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dw_dov_L1_claude.log

echo "### D_W + D_ov spectra  L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  config=${CONFIG}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

# ---- D_W ----
BIN_DW="eig_dw_spectrum_L${LREF}_nt${NT}.o"
echo "### BUILD (-DSPECTRUM_DW) -> $BIN_DW ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DSPECTRUM_DW -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DW" 2>>"$LOG" \
  || { echo "### D_W BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }
echo "### RUN D_W geev  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DW" 0.0 0.0 2>&1 | tee -a "$LOG" \
  || { echo "### D_W RUN FAILED ###" | tee -a "$LOG"; exit 1; }

# ---- D_ov ----
BIN_DOV="eig_dov_spectrum_L${LREF}_nt${NT}.o"
echo "### BUILD (-DSPECTRUM_DOV) -> $BIN_DOV ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DSPECTRUM_DOV -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DOV" 2>>"$LOG" \
  || { echo "### D_ov BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }
echo "### RUN D_ov geev  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DOV" 0.0 0.0 2>&1 | tee -a "$LOG" \
  || { echo "### D_ov RUN FAILED ###" | tee -a "$LOG"; exit 1; }

# ---- side-by-side plot ----
echo "### PLOT D_W vs D_ov  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 eig_dw_dov_compare_claude.py \
  "eig_dw_spectrum_L${LREF}_nt${NT}_claude.dat" \
  "eig_dov_spectrum_L${LREF}_nt${NT}_claude.dat" 2>&1 | tee -a "$LOG"

echo "### D_W + D_ov done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
