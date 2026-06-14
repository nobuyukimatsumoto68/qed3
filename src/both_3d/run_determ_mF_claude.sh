#!/usr/bin/env bash
# DETERMINISTIC ground truth at m_F = 0.1 (real, flavor-breaking), free U=1, L=1.  GPU 0.
# Pipeline (exact1 single-insertion for exact-K; SUMMED local; condensates).  Steps 1-4 READ the dense
# propagator from step 0, so step 0 MUST run first.
#   0) jj_propagator_deter        -> data_free_vmRe0.1.../prop_deter_L1/Dinv.0.h5  (dense D_m^{-1}+Dov)
#   1) jj_exact_diag_deter_free   --ins 0 -> corr_deter_exact1_L1         (vector exact-K tp/sp; Vmm=conj)
#   2) jj_exact_axial_deter_free  --ins 0 -> corr_deter_exact1_axial_L1   (axial exact-K, D_m legs)
#   3) jj_local_deter             (summed) -> corr_deter_local_L1         (local vector s1,s2,s3)
#   4) jj_local_axial_deter       (summed) -> corr_deter_local_axial_L1   (local axial s1,s2,s3)
#   5) condensate_deter_free      -> condensate_deter_L1                  (PS/FS)
# m_F is FULLY correct here (no tilde-D; Vmm=conj(Vpp), axial uses D_m -- matches the dilute's m_F path).
# Distinct binaries (*_detF.o) -> no clash with prod / other validation builds.  L=1 via -DN_REFINE_CLI=1.
# Run BEFORE run_jj_dilute_valid_mF (the notebook compares the stochastic run against these).  Fast (~minutes).
# Read back: determ_mF_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=determ_mF_claude.log
MRE=0.1; MIM=0.0

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
DEF="-DN_REFINE_CLI=1"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

: > "$LOG"
build(){   # $1 = source, $2 = output binary
  echo "### compile $1 -> $2  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  $NVCC $NVCCFLAGS $DEF $INCLUDES $LDFLAGS "$1" -o "$2" 2>&1 | tee -a "$LOG"
  if [ "${PIPESTATUS[0]}" -ne 0 ]; then echo "### BUILD FAILED: $1 ###" | tee -a "$LOG"; exit 1; fi
}
runstep(){ # $1 = label, rest = command
  local label="$1"; shift
  echo "### run $label  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  "$@" 2>&1 | tee -a "$LOG"
  if [ "${PIPESTATUS[0]}" -ne 0 ]; then echo "### STEP FAILED: $label ###" | tee -a "$LOG"; exit 1; fi
}

build jj_propagator_deter_claude.cu       jj_prop_detF.o
build jj_exact_diag_deter_free_claude.cu  jj_exdiag_detF.o
build jj_exact_axial_deter_free_claude.cu jj_exaxial_detF.o
build jj_local_deter_claude.cu            jj_local_detF.o
build jj_local_axial_deter_claude.cu      jj_localaxial_detF.o
build condensate_deter_free_claude.cu     condensate_detF.o

runstep "0 propagator (dense D_m^{-1}+Dov)" ./jj_prop_detF.o       --mass-re "$MRE" --mass-im "$MIM"
runstep "1 exact1 vector"                   ./jj_exdiag_detF.o     --mass-re "$MRE" --mass-im "$MIM" --ins 0 --n-t0 2
runstep "2 exact1 axial"                    ./jj_exaxial_detF.o    --mass-re "$MRE" --mass-im "$MIM" --ins 0 --n-t0 2
runstep "3 local vector (summed)"           ./jj_local_detF.o      --mass-re "$MRE" --mass-im "$MIM" --n-t0 2
runstep "4 local axial (summed)"            ./jj_localaxial_detF.o --mass-re "$MRE" --mass-im "$MIM" --n-t0 2
runstep "5 condensate"                      ./condensate_detF.o    --mass-re "$MRE" --mass-im "$MIM"

echo "### determ m_F=0.1 complete -> data_free_vmRe0.100000vmIm0.000000/{prop_deter,corr_deter_*,condensate_deter}_L1 ###" | tee -a "$LOG"
