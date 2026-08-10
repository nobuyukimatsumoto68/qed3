#!/usr/bin/env bash
# tmp_build_L3_claude.sh -- PRE-FLIGHT build of the NEW L=3 conn+disc measurement binaries, to run
# AFTER the machine/livepatch update and BEFORE the long conn->disc resume.  Verifies -DN_REFINE_CLI=3
# compiles and the n3 geometry is in place, so a build failure does NOT abort the full resume
# (run_conn_disc_ext_claude.sh does `exit 1` on any BUILD FAILED).  CPU compile only, no GPU, no rm.
#
# Run:  bash tmp_build_L3_claude.sh   (reads back tmp_build_L3_claude.log)
# If the machine's NVIDIA driver/CUDA changed in the update, also rebuild L1/L2/L4 on the resume with
#   FORCE_BUILD=1 nohup bash run_conn_disc_ext_claude.sh >> run_conn_disc_ext_claude.log 2>&1 &
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
LOG=tmp_build_L3_claude.log
: > "$LOG"

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

{
  echo "======== L3 pre-flight build $(date) ========"

  echo "---- check n3 geometry present ----"
  miss=0
  for f in pts nns links face alpha omega dual_links nns_dual face_dual alpha_dual omega_dual dualtriangleareas pts_dual
  do
    if [ -f "../../geometry/data/${f}_n3.dat" ]
    then
      echo "  OK ${f}_n3.dat"
    else
      echo "  MISSING ${f}_n3.dat"
      miss=1
    fi
  done
  [ "$miss" -eq 0 ] || { echo "### ABORT: n3 geometry incomplete -- (re)run geometry/geom.out 3 ###"; exit 1; }

  echo "---- compile conn L3 ----"
  $NVCC $NVCCBASE -DN_REFINE_CLI=3 $INCLUDES $LDFLAGS jj_local_ylm_scalar_conn_stoch_claude.cu -o jj_local_ylm_scalar_conn_stoch_L3.o \
    || { echo "### conn L3 BUILD FAILED ###"; exit 1; }
  echo "  conn L3 OK -> jj_local_ylm_scalar_conn_stoch_L3.o"

  echo "---- compile disc L3 ----"
  $NVCC $NVCCBASE -DN_REFINE_CLI=3 $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o jj_local_ylm_disc_stoch_L3.o \
    || { echo "### disc L3 BUILD FAILED ###"; exit 1; }
  echo "  disc L3 OK -> jj_local_ylm_disc_stoch_L3.o"

  echo "======== L3 pre-flight build DONE (both binaries OK) $(date) ========"
} 2>&1 | tee -a "$LOG"
