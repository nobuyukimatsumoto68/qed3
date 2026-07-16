#!/usr/bin/env bash
# tmp_eig_dov_spectrum_L1_scan_claude.sh
# _claude: COUPLING SCAN of the dense complex D_ov spectrum at L1 (N=3072), gsq = 2, 4, 8, each on a
# thermalized config.  Builds eig_wmass_val_claude.cu once with -DSPECTRUM_DOV; runs geev per gsq with an
# OUT_TAG so the three .dat files coexist; then a 3-panel comparison plot.  Shows whether the empty near-zero
# (left) arc of the GW circle fills in as the coupling weakens.  GPU 1.  Reads back:
# eig_dov_spectrum_L1_scan_claude.log + eig_dov_spectrum_compare_claude.png .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128

# gsq | thermalized config k
GSQS=(2 4 8)
KS=(4000 300 2105)

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dov_spectrum_L1_scan_claude.log

echo "### D_ov spectrum COUPLING SCAN  L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  gsq=2,4,8  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN="eig_dov_spectrum_L${LREF}_nt${NT}.o"
echo "### BUILD (-DSPECTRUM_DOV) -> $BIN ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DSPECTRUM_DOV -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN" 2>>"$LOG" \
  || { echo "### BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }

for i in "${!GSQS[@]}"; do
  g="${GSQS[$i]}"
  k="${KS[$i]}"
  cfg="Nf2_gsq${g}.000000at0.200000nu01.000000nt128L1/ckpoint_lat.${k}"
  echo "### RUN gsq=${g}  config k=${k}  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  CONFIG_LAT="$cfg" OUT_TAG="_gsq${g}" CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" 0.0 0.0 2>&1 | tee -a "$LOG" \
    || { echo "### gsq=${g} RUN FAILED ###" | tee -a "$LOG"; exit 1; }
done

echo "### COMPARISON PLOT  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 eig_dov_spectrum_compare_claude.py 2>&1 | tee -a "$LOG"

echo "### D_ov coupling scan done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
