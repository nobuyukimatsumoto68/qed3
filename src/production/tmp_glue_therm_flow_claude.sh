#!/usr/bin/env bash
# tmp_glue_therm_flow_claude.sh -- build the PRODUCTION h5 flow driver (L1/L2/L4) and smoke-test:
# (1) run L1 g0.5 at stride 100 (~20 configs), (2) RERUN to verify resume (must flow nothing),
# (3) print the h5 shapes. Protocol: eps=0.01, save every step, tmax=2.0 -> tlist 201 points.
# Run:  bash tmp_glue_therm_flow_claude.sh       (tee'd to tmp_glue_therm_flow_claude.log)
# Writes: glue_therm_flow_L{1,2,4}_claude.o and
#   data_Nf2_gsq0.500000..._hb1.000000/therm_flow_claude.h5  (smoke; kept -- production reuses it)
# No rm anywhere in this script. Configs read-only.
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
export OMP_NUM_THREADS=4

LOG=tmp_glue_therm_flow_claude.log

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lm'
SRC=glue_therm_flow_claude.cu

{
  echo "##################################################################"
  echo "### build glue_therm_flow (L1/L2/L4)  $(date)"
  echo "##################################################################"
} | tee "$LOG"

for L in 1 2 4
do
  BIN=glue_therm_flow_L${L}_claude.o
  echo "### building $BIN ###" | tee -a "$LOG"
  $NVCC $NVCCBASE $INCLUDES $LDFLAGS -DN_REFINE_CLI=$L -o "$BIN" "$SRC" 2>&1 | tee -a "$LOG"
  rc=${PIPESTATUS[0]}
  if [ "$rc" -ne 0 ]
  then
    echo "### ERROR: build failed for L=$L (exit $rc) -- aborting ###" | tee -a "$LOG"
    exit "$rc"
  fi
done

ENS=Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000

{
  echo "##################################################################"
  echo "### smoke 1/2: L1 g0.5 stride 100 (fresh or resume)  $(date)"
  echo "##################################################################"
} | tee -a "$LOG"
./glue_therm_flow_L1_claude.o 0.5 2 1.0 100000 1 100 "$ENS"/ 2.0 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: smoke run failed (exit $rc) ###" | tee -a "$LOG"
  exit "$rc"
fi

{
  echo "##################################################################"
  echo "### smoke 2/2: RERUN -- must resume with 0 new flows  $(date)"
  echo "##################################################################"
} | tee -a "$LOG"
./glue_therm_flow_L1_claude.o 0.5 2 1.0 100000 1 100 "$ENS"/ 2.0 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]
then
  echo "### ERROR: resume rerun failed (exit $rc) ###" | tee -a "$LOG"
  exit "$rc"
fi

{
  echo "### h5 content check ###"
  python3 - <<PY
import h5py
f = h5py.File("data_$ENS/therm_flow_claude.h5", "r")
for name in ["tlist", "klist", "E", "eps", "save_every", "tmax"]:
    d = f[name]
    print(name, d.shape, d[()].ravel()[:3], "...")
PY
  echo "### done  $(date) ###"
} | tee -a "$LOG"
