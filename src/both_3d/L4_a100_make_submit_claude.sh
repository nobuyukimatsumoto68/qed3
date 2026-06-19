#!/bin/bash
# L4_a100_make_submit_claude.sh -- compile the L=4 free-propagator binaries (sm_80, A100) and sbatch the job.
#   builds : L4_a100_prop_4.o      (stage 1, A100 dense overlap propagator -- needs cuSOLVER)
#            jj_local_deter_4.o    (stages 2 & 3, host local-op contraction: lattice P and continuum P)
#   then   : sbatch [--dependency=afterany:<DEP_JOBID>] L4_a100_batch_claude.sh
# Compile flags mirror Makefile.fnal (cluster HDF5/HighFive/qfe_mod/gsl include + lib paths) PLUS
# -lcusolver -lcusparse, which the propagator needs (cusolverDnZgetrf/Zgetrs/Xgeev) and Makefile.fnal omits.
# L (= N_REFINE) is COMPILE-TIME: -DN_REFINE_CLI=4 ; Nt defaults to 128 (override with -DNT_CLI=...).
# Run from src/both_3d/:  bash L4_a100_make_submit_claude.sh [DEP_JOBID]
#   DEP_JOBID (optional): queue the batch with --dependency=afterany:DEP_JOBID so it starts only AFTER that
#   job ENDS (any exit, incl. TIMEOUT) -- use to submit the follow-on NOW behind a still-running stage-2 job
#   (e.g. 1286235).  afterany (NOT afterok): a TIMEOUT is a non-zero exit, and afterok would then be marked
#   DependencyNeverSatisfied and never run.
# Compiles to <bin>.new then atomically mv's into place, so a binary the running job is currently executing
# is not clobbered mid-run ('text file busy').
# set -e: stop immediately on any failure (incl. a failed compile) so the job is NOT submitted on a bad build.
set -eu

cd /project/qed3/qed3/src/both_3d || { echo "cd to src/both_3d FAILED"; exit 1; }
source /lustre2/qed3/env.sh

L=4
DEP_JOBID="${1:-}"   # optional: job id to wait for (afterany); empty => submit immediately

NVCC=nvcc
# CUDA_HOME/CUDA_MATH/CUDA_COMP may be unset after env.sh; with set -u, mirror Makefile.fnal's empty-expansion
# via ${VAR:+ -L.../lib...} so an unset var contributes nothing (the modules already put the libs on the path).
NVCCFLAGS="-arch=sm_80 -g -O3 -std=c++20 ${CUDA_HOME:+-L${CUDA_HOME}/lib64 -L${CUDA_HOME}/lib64/stubs} ${CUDA_MATH:+-L${CUDA_MATH}/lib64} ${CUDA_COMP:+-L${CUDA_COMP}/lib} -lcudart -lcuda -lnvidia-ml -lcublas -lcusolver -lcusparse -lcufft -ldl -lgomp -Xcompiler -fopenmp"
INCLUDES='-I/project/qed3/qed3/src/both_3d/includes/ -I/project/qed3/qed3/qfe_mod/include/ -I/project/qed3/highfive/include/ -I/srv/software/el8/x86_64/eb/HDF5/1.14.2-GCC-12.3.0-serial/include/ -I/project/qed3/gsl/include/'
LDFLAGS='-L/srv/software/el8/x86_64/eb/HDF5/1.14.2-GCC-12.3.0-serial/lib/ -L/project/qed3/gsl/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "===== compile stage 1: L4_a100_prop (N_REFINE_CLI=${L}, sm_80, +cusolver) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} L4_a100_prop_claude.cu   -o L4_a100_prop_${L}.o.new \
  || { echo "BUILD FAILED (prop)";  exit 1; }
mv -f L4_a100_prop_${L}.o.new L4_a100_prop_${L}.o   # atomic; safe even if a running job holds the old inode

echo "===== compile stages 2&3: jj_local_deter (N_REFINE_CLI=${L}, sm_80) ====="
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_local_deter_claude.cu -o jj_local_deter_${L}.o.new \
  || { echo "BUILD FAILED (local)"; exit 1; }
mv -f jj_local_deter_${L}.o.new jj_local_deter_${L}.o   # atomic; running 1286235 keeps its old inode

if [ -n "${DEP_JOBID}" ]; then
  echo "===== submit L4_a100_batch_claude.sh  (dependency afterany:${DEP_JOBID}) ====="
  sbatch --dependency=afterany:${DEP_JOBID} L4_a100_batch_claude.sh
else
  echo "===== submit L4_a100_batch_claude.sh ====="
  sbatch L4_a100_batch_claude.sh
fi
