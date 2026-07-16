#!/usr/bin/env bash
# tmp_eig_dw_dov_L1_scan_claude.sh
# _claude: D_W vs D_ov side-by-side (unitary-projector picture) at L1 for gsq = 2 and 4 (gsq8 already done),
# each on a thermalized config.  Shows how the Wilson eigenvalue RING moves with coupling and whether it
# approaches the zero-crossing direction (which would create near-zero overlap modes).  Two builds of
# eig_wmass_val_claude.cu (-DSPECTRUM_DW, -DSPECTRUM_DOV), reused across couplings via CONFIG_LAT/OUT_TAG;
# one 2-panel figure per gsq.  N=3072 (~1 GB).  GPU 1.  Reads back: eig_dw_dov_L1_scan_claude.log +
# eig_dw_dov_compare_L1_gsq{2,4}_claude.png .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
LREF=1
NT=128

# gsq | thermalized config k
GSQS=(2 4)
KS=(4000 300)

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dw_dov_L1_scan_claude.log

echo "### D_W vs D_ov COUPLING SCAN  L=${LREF} Nt=${NT} (N=$((2*(10*LREF*LREF+2)*NT)))  gsq=2,4  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

BIN_DW="eig_dw_spectrum_L${LREF}_nt${NT}.o"
BIN_DOV="eig_dov_spectrum_L${LREF}_nt${NT}.o"
echo "### BUILD D_W + D_ov ###" | tee -a "$LOG"
$NVCC $NVCCBASE -DSPECTRUM_DW  -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DW"  2>>"$LOG" || { echo "### D_W BUILD FAILED ###" | tee -a "$LOG"; exit 1; }
$NVCC $NVCCBASE -DSPECTRUM_DOV -DLREF=${LREF} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_wmass_val_claude.cu -o "$BIN_DOV" 2>>"$LOG" || { echo "### D_ov BUILD FAILED ###" | tee -a "$LOG"; exit 1; }

for i in "${!GSQS[@]}"; do
  g="${GSQS[$i]}"
  k="${KS[$i]}"
  cfg="Nf2_gsq${g}.000000at0.200000nu01.000000nt128L1/ckpoint_lat.${k}"
  tag="_gsq${g}"
  echo "### gsq=${g}  config k=${k}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  CONFIG_LAT="$cfg" OUT_TAG="$tag" CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DW"  0.0 0.0 2>&1 | tee -a "$LOG" || { echo "### gsq=${g} D_W RUN FAILED ###" | tee -a "$LOG"; exit 1; }
  CONFIG_LAT="$cfg" OUT_TAG="$tag" CUDA_VISIBLE_DEVICES=$GPU ./"$BIN_DOV" 0.0 0.0 2>&1 | tee -a "$LOG" || { echo "### gsq=${g} D_ov RUN FAILED ###" | tee -a "$LOG"; exit 1; }
  python3 eig_dw_dov_compare_claude.py \
    "eig_dw_spectrum_L${LREF}${tag}_nt${NT}_claude.dat" \
    "eig_dov_spectrum_L${LREF}${tag}_nt${NT}_claude.dat" \
    "eig_dw_dov_compare_L${LREF}${tag}_claude.png" \
    "(L1, gsq=${g})" 2>&1 | tee -a "$LOG"
done

echo "### D_W vs D_ov scan done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
