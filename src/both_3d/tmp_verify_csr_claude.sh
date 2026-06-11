#!/usr/bin/env bash
# One-time verification of the O(len) CSR bucketing (sparse_dirac_claude.h) against the ORIGINAL
# O(N*len) scan.  Builds jj_disp_deter at L=1 with -DCSR_VERIFY: DWDevice::initialize() rebuilds the
# old scan into check arrays and assert()s it equals the bucketing, printing
#     # CSR_VERIFY: O(len) bucketing == O(N*len) scan  (N=..., nnz=...)
# on success (the assert aborts on any mismatch).  The check + the bucketing both run at overlap
# construction (startup), BEFORE any propagator load or output.
#
# CHEAPEST + NON-DESTRUCTIVE:
#   - L=1 (N=3072): both the bucketing and the slow old scan are instant.
#   - jj_disp_deter is the lightest consumer (no dense-K build).
#   - --out-tag verify => output goes ONLY to the throwaway dir corr_deter_disp_verify_L1/ ;
#     no existing corr_deter_* data is read-for-overwrite or touched.  (If that dir already exists and
#     is complete, the run SKIPS the compute but STILL prints CSR_VERIFY at construction -- fine.)
#   - Requires the L=1 lattice propagator data_free_*/prop_deter_L1/Dinv.0.h5 to exist (it does from
#     earlier disp/loc runs); only read, never written.
#
# NOTE: asserts must be enabled -- the project's NVCCFLAGS do NOT set -DNDEBUG, so they are.  Seeing the
# "CSR_VERIFY: ... ==" line = all asserts passed.  A separate verify binary name is used so the real
# jj_disp_deter_L1.o is left intact.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp -DCSR_VERIFY"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=csr_verify_L1_claude.log
echo "### compile jj_disp_deter L=1 with -DCSR_VERIFY ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=1 jj_disp_deter_claude.cu -o jj_disp_deter_verify_L1.o \
  2>&1 | tee -a "$LOG" || { echo "BUILD FAILED"; exit 1; }
echo "### run L=1 (--out-tag verify => only corr_deter_disp_verify_L1/ is written) ###" | tee -a "$LOG"
./jj_disp_deter_verify_L1.o --out-tag verify --n-t0 2 \
  2>&1 | tee -a "$LOG" || { echo "RUN FAILED"; exit 1; }
echo "### done -- success = the 'CSR_VERIFY: O(len) bucketing == O(N*len) scan' line above ###" | tee -a "$LOG"
