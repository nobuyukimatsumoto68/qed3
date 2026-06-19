#!/usr/bin/env bash
# _claude: build + run the measure-weighted diagonal-mass L=1/L=2 check.
# Validates: (i) operator obsolete-vs-production at L=1 (mult/adj/DHD + _ms), (ii) HMC force
# obsolete-vs-production at L=1, (iii) force-vs-finite-difference at L=1 and L=2, for BOTH the
# default reference grad and the active block grad_l4 (-DGRAD_L4). Nontrivial gaussian gauge.
# Reads back via test_diag_mass_claude.log. No rm anywhere (per project rule).
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=test_diag_mass_claude.log
GPU=0                       # MPS is up on either GPU; set to 1 to use the other one
export CUDA_VISIBLE_DEVICES=$GPU

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_diag_mass_l1_claude.cu

run_phase () {
  local tag="$1"
  local lref="$2"
  local extra="$3"
  local bin="test_diag_${tag}.x"
  echo ""
  echo "==================== BUILD ${tag}  (LREF=${lref} ${extra}) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=${lref} ${extra} $SRC $LDFLAGS -o "$bin" ; then
    echo "!!! BUILD FAILED: ${tag} -- skipping run"
    return
  fi
  echo "-------------------- RUN ${tag} (GPU ${GPU}) --------------------"
  if ./"$bin" ; then
    echo ">>> ${tag}: exit 0 (ALL PASS)"
  else
    echo ">>> ${tag}: NONZERO exit (some check FAILED -- inspect above)"
  fi
}

{
  echo "######## diagonal-mass check  $(date)  GPU=${GPU} ########"
  # All four grad variants at L=1 + L=2 (default / l1 / l2 / l4); default + l4 at L=4.
  # L=1: obsolete-vs-production force at machine precision (catches a wrong mass_coeff/conj
  # in ANY grad variant, even though m_L is uniform there). L=2,4: site-varying m_L via force-vs-FD.
  run_phase "L1_default" 1 ""
  run_phase "L1_gradl1"  1 "-DGRAD_L1"
  run_phase "L1_gradl2"  1 "-DGRAD_L2"
  run_phase "L1_gradl4"  1 "-DGRAD_L4"
  run_phase "L2_default" 2 ""
  run_phase "L2_gradl1"  2 "-DGRAD_L1"
  run_phase "L2_gradl2"  2 "-DGRAD_L2"
  run_phase "L2_gradl4"  2 "-DGRAD_L4"
  run_phase "L4_default" 4 ""
  run_phase "L4_gradl4"  4 "-DGRAD_L4"
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
