#!/usr/bin/env bash
# tmp_eig_dw_dov_L1_free_claude.sh
# _claude: D_W vs D_ov side-by-side at L1 in the FREE field (U=0, cold gauge) -- the clean reference with no
# gauge broadening, so the doubler branches are sharp and the domain wall M=1 sits at its ideal location.
# Same two builds of eig_wmass_val_claude.cu (-DSPECTRUM_DW, -DSPECTRUM_DOV); CONFIG_LAT is NOT set, so the
# driver uses the cold gauge.  N=3072.  GPU 1.  Reads back: eig_dw_dov_L1_free_claude.log +
# eig_dw_dov_compare_L1_free_claude.png .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128
unset CONFIG_LAT          # free field (cold gauge)
export OUT_TAG="_free"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dw_dov_L1_free_claude.log

echo "### D_W vs D_ov  L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  FREE FIELD  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN_DW="eig_dw_spectrum_L${LREF}_nt${NT}.o"
BIN_DOV="eig_dov_spectrum_L${LREF}_nt${NT}.o"
echo "### BUILD D_W + D_ov ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DSPECTRUM_DW  -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DW"  2>>"$LOG" || { echo "### D_W BUILD FAILED ###" | tee -a "$LOG"; exit 1; }
$NVCC $NVCCBASE -DSPECTRUM_DOV -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DOV" 2>>"$LOG" || { echo "### D_ov BUILD FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### RUN D_W (free)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DW"  0.0 0.0 2>&1 | tee -a "$LOG" || { echo "### D_W RUN FAILED ###" | tee -a "$LOG"; exit 1; }
echo "### RUN D_ov (free)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DOV" 0.0 0.0 2>&1 | tee -a "$LOG" || { echo "### D_ov RUN FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### PLOT  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 eig_dw_dov_compare_claude.py \
  "eig_dw_spectrum_L${LREF}_free_nt${NT}_claude.dat" \
  "eig_dov_spectrum_L${LREF}_free_nt${NT}_claude.dat" \
  "eig_dw_dov_compare_L${LREF}_free_claude.png" \
  "(L1, free field)" 2>&1 | tee -a "$LOG"

echo "### free field done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
