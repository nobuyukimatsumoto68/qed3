#!/bin/bash -l
set -u

# ============================================================================
# Handoff: build + run the SHAPE-BASIS glueball measurements on the massless
# L=1 gsq8 SEA ensembles (Nf2/4/6). Run this one script; read back the logs.
#
#   glue_f2_shapes_claude.cu    -> F^2 / 0++ (squared), ell 0-3, dir data_..._shapes_f2/
#   glue2_msm_shapes_claude.cu  -> linear F_12,          ell 1-3, dir data_..._shapes_msm/
#
# Basis = 5 icosahedral-orbit shapes {triangle, rect, twisted-rect, figure-8,
# twisted-figure-8}; single flow tmax=1.0/100. Config selection kmin=100 stride=1.
# 6 single-threaded processes (2 drivers x 3 Nf) -> 6 cores. Resume-safe (h5 complete flag).
# NO rm anywhere.
#
# *** IMPORTANT ***  The Nf2 data_..._shapes_{f2,msm}/ dirs currently hold 11 STALE
# smoke files (older operator count). Because resume only checks the complete flag,
# NOT nops, a dir must be internally consistent. rm the stale smoke dirs BY HAND first:
#   rm -r data_Nf2_gsq8.000000at0.200000nu01.000000nt128L1_shapes_f2
#   rm -r data_Nf2_gsq8.000000at0.200000nu01.000000nt128L1_shapes_msm
# (Nf4/Nf6 dirs are fresh; no action needed there.)
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

# module load gcc/13.2.0
# module load cuda/12.6

export OMP_NUM_THREADS=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

GSQ=8.0
NU0=1.0
KMAXRUN=100000
KMIN=100
STRIDE=1

BLOG=glue_shapes_L1_build_claude.log

echo "================ BUILD glue2_msm_shapes_claude.cu (linear) ================" | tee "$BLOG"
date | tee -a "$BLOG"
"$NVCC" glue2_msm_shapes_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_claude.o 2>&1 | tee -a "$BLOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "linear BUILD FAILED rc=$rc" | tee -a "$BLOG"; exit "$rc"; fi

echo "================ BUILD glue_f2_shapes_claude.cu (F^2) ================" | tee -a "$BLOG"
date | tee -a "$BLOG"
"$NVCC" glue_f2_shapes_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_claude.o 2>&1 | tee -a "$BLOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "F^2 BUILD FAILED rc=$rc" | tee -a "$BLOG"; exit "$rc"; fi
echo "BUILD OK (both)" | tee -a "$BLOG"

echo "================ LAUNCH 6 measurement processes ================" | tee -a "$BLOG"
for NF in 2 4 6; do
  ./glue_f2_shapes_claude.o   "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > glue_f2_shapes_Nf${NF}_run_claude.log  2>&1 &
  ./glue2_msm_shapes_claude.o "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > glue_msm_shapes_Nf${NF}_run_claude.log 2>&1 &
done
echo "launched 6 procs; waiting..." | tee -a "$BLOG"
wait

echo "================ ALL DONE ================" | tee -a "$BLOG"
date | tee -a "$BLOG"
for tag in shapes_f2 shapes_msm; do
  for nf in 2 4 6; do
    d="data_Nf${nf}_gsq8.000000at0.200000nu01.000000nt128L1_${tag}"
    echo "  $tag Nf$nf : $(ls "$d"/F_corr.*.h5 2>/dev/null | wc -l) F_corr.h5 files" | tee -a "$BLOG"
  done
done
