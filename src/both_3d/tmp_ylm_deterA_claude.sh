#!/usr/bin/env bash
# Chunk A: deterministic per-m Ylm tower GROUND TRUTH (jj_local_deter_claude.cu, generalized to
# 3 Pauli channels s1/s2/s3 and l<=3, per-m output).  Free field, L=1, LATTICE overlap propagator
# (reads the EXISTING data_free_.../prop_deter_L1/Dinv.0.h5 -- no LU rebuild), --n-t0 2.
#
# Uses --out-tag ylmA so it writes a FRESH dir  data_free_vmRe0.000000vmIm0.000000/corr_deter_local_ylmA_L1/
# -- it does NOT touch the existing corr_deter_local_L1 (complete-gated old keys), so NO rm needed.
# Distinct binary jj_local_deter_ylmA_L1.o (no ETXTBSY with other live runs).
#
# New ylm keys (per-m, 3 Pauli):  h0/t0_b/ylm/s{a}/l{0..3}/m{-l..l}/{Vpp[,Vmm]}   (a=1,2,3).
# The s1/s2/s3 local + disc keys are UNCHANGED from before.  Read back: jj_ylm_deterA_claude.log .
# CHECK (Chunk C notebook): tp=s3 rates ->(2,3,4), G22 e^{3t}/G11 e^{2t}->2.4, g0->0, sp/tp->-1.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
L=1
BIN=jj_local_deter_ylmA_L${L}.o
LOG=jj_ylm_deterA_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

echo "### compile -> $BIN  L=${L}  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} jj_local_deter_claude.cu -o "$BIN" 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run free L=${L} (lattice prop prop_deter_L1), --n-t0 2 --out-tag ylmA  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./"$BIN" --n-t0 2 --out-tag ylmA 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### RUN FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### done -> data_free_vmRe0.000000vmIm0.000000/corr_deter_local_ylmA_L1/corr.0.h5  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
