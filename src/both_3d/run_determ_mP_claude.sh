#!/usr/bin/env bash
# DETERMINISTIC ground truth at m_P = 0.1 i (imaginary, parity-breaking), free U=1, L=1.  GPU 1.
#
# tilde-D support ADDED (2026-06-13): jj_propagator_deter saves Dtil_inv = tilde D_{m_P}^{-1}; the exact-K
# deter now uses it for the parity-asymmetric "-" channels (Eqs. 3.65/3.67):
#   - vector Vpp  : K . D_m^{-1}                                  (Eq. 3.63)      VALID
#   - vector Vmm  : K^dag . tilde D_{m_P}^{-dag}  (jj_exact_diag) (Eq. 3.65)      VALID (NEW)
#   - axial  Apm  : t0-leg tilde D^{-dag}, sink D_m^{-1} (jj_exact_axial, Eq. 3.67) VALID (NEW)
#   - condensate  : condensate_deter_free (tilde + (1+m_P)^{-1})                  VALID
# LOCAL channels:
#   - local vector Vpp : D_m  VALID ; local vector Vmm : SKIPPED at parity in BOTH (nothing to compare).
#   - local axial      : tilde-D t0-leg + (1+m_P)^{-1} (Eq. 3.61) in BOTH the dilute AND jj_local_axial_deter
#       (the bare current has no conserved-kernel correction, so the (1+m_P)^{-1} is NOT cancelled, unlike
#       the exact-K axial Eq. 3.67).  -> VALID (fixed 2026-06-13).
#
# Distinct binaries (*_detP.o) -> no clash.  L=1 via -DN_REFINE_CLI=1.  Run BEFORE run_jj_dilute_valid_mP.
# Read back: determ_mP_claude.log
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=determ_mP_claude.log
MRE=0.0; MIM=0.1

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

build jj_propagator_deter_claude.cu       jj_prop_detP.o
build jj_exact_diag_deter_free_claude.cu  jj_exdiag_detP.o
build jj_exact_axial_deter_free_claude.cu jj_exaxial_detP.o
build jj_local_deter_claude.cu            jj_local_detP.o
build jj_local_axial_deter_claude.cu      jj_localaxial_detP.o
build condensate_deter_free_claude.cu     condensate_detP.o

runstep "0 propagator (D_m^{-1}+Dtil_inv+Dov)"     ./jj_prop_detP.o       --mass-re "$MRE" --mass-im "$MIM"
runstep "1 exact1 vector (Vpp + Vmm[tilde])"        ./jj_exdiag_detP.o     --mass-re "$MRE" --mass-im "$MIM" --ins 0 --n-t0 2
runstep "2 exact1 axial (tilde t0-leg)"             ./jj_exaxial_detP.o    --mass-re "$MRE" --mass-im "$MIM" --ins 0 --n-t0 2
runstep "3 local vector (Vpp; Vmm skipped, matches dilute)" ./jj_local_detP.o --mass-re "$MRE" --mass-im "$MIM" --n-t0 2
runstep "4 local axial (tilde t0-leg + (1+m_P)^{-1}; Eq. 3.61)" ./jj_localaxial_detP.o --mass-re "$MRE" --mass-im "$MIM" --n-t0 2
runstep "5 condensate (full m_P support)"           ./condensate_detP.o    --mass-re "$MRE" --mass-im "$MIM"

echo "### determ m_P=0.1i complete.  VALID: vector Vpp+Vmm + axial (exact-K, tilde), local Vpp + local axial" | tee -a "$LOG"
echo "### (tilde t0-leg + (1+m_P)^{-1}), condensate.  Local vector Vmm skipped both sides (nothing to compare). ###" | tee -a "$LOG"
