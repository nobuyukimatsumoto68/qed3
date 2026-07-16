#!/usr/bin/env bash
# tmp_eig_dov_spectrum_L1_claude.sh
# _claude: DENSE complex spectrum of the overlap operator D_ov (NOT D_ov^dag D_ov) at L1 on a THERMALIZED
# config -> the eigenvalues should sit on the GW circle |z-1|=1.  Cheap at L1 (N=3072).  Builds
# eig_wmass_val_claude.cu with -DSPECTRUM_DOV (binds mult, writes Re/Im/|z|), runs geev, then plots the
# complex plane.  Single process on GPU 1.  Reads back: eig_dov_spectrum_L1_claude.log + the .dat + PNG.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128
CONFIG="Nf2_gsq8.000000at0.200000nu01.000000nt128L1/ckpoint_lat.2105"   # thermalized (traj 2105/3426); smallest
                                                                       # Wilson lambda_min among thermalized L1
export CONFIG_LAT="$CONFIG"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dov_spectrum_L1_claude.log

echo "### D_ov COMPLEX spectrum  L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  config=${CONFIG}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN="eig_dov_spectrum_L${LREF}_nt${NT}.o"
echo "### BUILD (-DSPECTRUM_DOV) -> $BIN ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DSPECTRUM_DOV -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN" 2>>"$LOG" \
  || { echo "### BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

echo "### RUN geev (mass=0)  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" 0.0 0.0 2>&1 | tee -a "$LOG" \
  || { echo "### RUN FAILED ###" | tee -a "$LOG"; exit 1; }

echo "### PLOT complex plane  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 eig_dov_spectrum_plot_claude.py "eig_dov_spectrum_L${LREF}_nt${NT}_claude.dat" 2>&1 | tee -a "$LOG"

echo "### D_ov spectrum done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
