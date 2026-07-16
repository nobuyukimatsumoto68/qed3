#!/usr/bin/env bash
# tmp_eig_dw_dov_L1_gsq1_claude.sh
# _claude: D_W vs D_ov side-by-side at L1, gsq=1 (weakest coupling) on a thermalized config.  Extends the
# gsq=2,4,8 scan -- expect the D_W ring's left edge to push further across the origin (Re<0) and the overlap
# arc to reach closer to z=0 (smaller |z|_min).  Two builds of eig_wmass_val_claude.cu (-DSPECTRUM_DW,
# -DSPECTRUM_DOV), one 2-panel figure.  N=3072 (~1 GB).  GPU 1.  Reads back: eig_dw_dov_L1_gsq1_claude.log +
# eig_dw_dov_compare_L1_gsq1_claude.png .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128
CONFIG="Nf2_gsq1.000000at0.200000nu01.000000nt128L1/ckpoint_lat.300"
export CONFIG_LAT="$CONFIG"
export OUT_TAG="_gsq1"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dw_dov_L1_gsq1_claude.log

echo "### D_W vs D_ov  L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  gsq=1  config=${CONFIG}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN_DW="eig_dw_spectrum_L${LREF}_nt${NT}.o"
BIN_DOV="eig_dov_spectrum_L${LREF}_nt${NT}.o"
echo "### BUILD D_W + D_ov ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DSPECTRUM_DW  -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DW"  2>>"$LOG" || { echo "### D_W BUILD FAILED ###" | tee -a "$LOG"; exit 1; }
$NVCC $NVCCBASE -DSPECTRUM_DOV -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DOV" 2>>"$LOG" || { echo "### D_ov BUILD FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### RUN D_W  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DW"  0.0 0.0 2>&1 | tee -a "$LOG" || { echo "### D_W RUN FAILED ###" | tee -a "$LOG"; exit 1; }
echo "### RUN D_ov  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DOV" 0.0 0.0 2>&1 | tee -a "$LOG" || { echo "### D_ov RUN FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### PLOT  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 eig_dw_dov_compare_claude.py \
  "eig_dw_spectrum_L${LREF}_gsq1_nt${NT}_claude.dat" \
  "eig_dov_spectrum_L${LREF}_gsq1_nt${NT}_claude.dat" \
  "eig_dw_dov_compare_L${LREF}_gsq1_claude.png" \
  "(L1, gsq=1)" 2>&1 | tee -a "$LOG"

echo "### gsq=1 done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
