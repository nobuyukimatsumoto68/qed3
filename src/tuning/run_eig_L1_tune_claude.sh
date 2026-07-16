#!/bin/bash
# _claude: settle the L=1 Wilson spectrum (lambda_min = smallest singular value of D_W, lambda_max) at
# the CORRECTED gsq range, using the fresh L1 tuning configs (gsq 0.5/1.0/1.5, 24 configs each). We have
# the gsq8 record (frozen_window_claude.md) but NOT these couplings -- the HMC only logs the IMPOSED
# frozen window (0.1,13), not the true spectrum. eig_wilson_lowmode_claude.cu measures it per config.
#
# The tuning dirs carry an "_hb1.000000" tag the eig dir-lookup doesn't try, so we symlink the bare name
# it expects -> the _hb dir. Output: wilson_lowmode_Nf2_gsq<g>_L1_claude.dat (cols: k lambda_min lambda_max
# ratio; header has the FREE-field U=0 reference). Light (N=3072, 24 configs) -> runs fast even under MPS.
#
# Run:  bash run_eig_L1_tune_claude.sh 2>&1 | tee run_eig_L1_tune_claude.log
set -u

SRCDIR=/mnt/barracuda22/qed3/qed3/src/tuning
cd "$SRCDIR" || exit 1
source ../../env.sh
export OMP_NUM_THREADS=4   # cap CPU threads (gpu0 already runs 2 L4mix jobs)

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=../both_3d/eig_wilson_lowmode_claude.cu   # measures spectrum UNFROZEN (compute_lambda_max)

# ---- build (LREF=1) -----------------------------------------------------------------------------
echo "===== build eig_L1.out (LREF=1)  [$(date +%F_%H:%M:%S)] ====="
$NVCC $NVCCFLAGS $INCLUDES -DLREF=1 "$SRC" $LDFLAGS -o eig_L1.out 2>&1 | tee build_eig_L1_claude.log
test -f eig_L1.out || { echo "BUILD FAILED"; exit 1; }

# ---- symlink the bare dir name the eig driver expects -> the _hb tuning dir ----------------------
AT=0.200000; NU0=1.000000; HB=_hb1.000000
for g in 0.500000 1.000000 1.500000; do
  base="Nf2_gsq${g}at${AT}nu0${NU0}mRe0.000000mIm0.000000nt128L1"
  hbdir="${base}${HB}"
  if [ -d "$hbdir" ] && [ ! -e "$base" ]; then ln -s "$hbdir" "$base"; echo "symlink $base -> $hbdir"; fi
done

# ---- scan each gsq (Nf=2, stride=1 = all 24 configs, kmax=0 = all) -------------------------------
for g in 0.5 1.0 1.5; do
  echo "===== eig scan L1 gsq$g  [$(date +%F_%H:%M:%S)] ====="
  CUDA_VISIBLE_DEVICES=0 ./eig_L1.out $g 2 1 0
done
echo "===== done -- see wilson_lowmode_Nf2_gsq{0.5,1.0,1.5}_L1_claude.dat  [$(date +%F_%H:%M:%S)] ====="
