#!/bin/bash -l
set -u

# ============================================================================
# Handoff: build + run BOTH glueball measurements on the massless L=1 gsq8 SEA
# ensembles (Nf2/4/6). Run this one script; read back the per-process logs.
#
#   glue_f2_claude.cu    -> F_{\mu\nu}^2 action density  (B^2/E^2, full Ylm tower incl l=0)
#   glue2_msm_claude.cu  -> linear F_12                   (Ylm tower l=1..3; l=0 dropped)
#
# Config selection: kmin=100 (skip thermalization), stride=1 (ALL configs), identical
#   config set for both operators. Resume-safe (skips configs whose F_corr already exists).
# Parallelism is ENSEMBLE-level: 6 single-threaded processes (2 drivers x 3 Nf) -> 6 cores.
#   No OMP over configs. Inner loops keep the existing NPARALLEL_GAUGE knob (=1).
# NO rm anywhere (project rule). A truncated F_corr from a previous killed run is NOT
#   auto-deleted; the PREFLIGHT below aborts and lists the file(s) for you to rm by hand.
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

# module load gcc/13.2.0   # uncomment on the cluster
# module load cuda/12.6

export OMP_NUM_THREADS=1   # each process is single-threaded; ensemble-level parallelism

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

GSQ=8.0
NU0=1.0
KMAXRUN=100000
KMIN=100
STRIDE=1
NT=128                  # expected F_corr row count (one row per dt)

BLOG=glue_massless_L1_build_claude.log

dir_of () { echo "data_Nf${1}_gsq8.000000at0.200000nu01.000000nt128L1_${2}"; }

# ---------------------------------------------------------------------------
echo "================ BUILD glue2_msm_claude.cu (F_12) ================" | tee "$BLOG"
date | tee -a "$BLOG"
"$NVCC" glue2_msm_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_claude.o 2>&1 | tee -a "$BLOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "F_12 BUILD FAILED rc=$rc" | tee -a "$BLOG"; exit "$rc"; fi

echo "================ BUILD glue_f2_claude.cu (F^2) ================" | tee -a "$BLOG"
date | tee -a "$BLOG"
"$NVCC" glue_f2_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_claude.o 2>&1 | tee -a "$BLOG"
rc=${PIPESTATUS[0]}
if [ "$rc" -ne 0 ]; then echo "F^2 BUILD FAILED rc=$rc" | tee -a "$BLOG"; exit "$rc"; fi
echo "BUILD OK (both)" | tee -a "$BLOG"

# ---------------------------------------------------------------------------
# Output is now per-config HDF5 (F_corr.<k>.h5) with a "complete" flag written last.
# Resume is robust: a killed/incomplete h5 lacks "complete", so it is re-measured
# (Truncate overwrites). No truncated-file preflight needed. The old text F_corr.k
# files (from the earlier text run) are ignored -- rm them by hand to reclaim disk.
echo "================ LAUNCH 6 measurement processes ================" | tee -a "$BLOG"
for NF in 2 4 6; do
  ./glue_f2_claude.o   "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > glue_f2_Nf${NF}_run_claude.log  2>&1 &
  ./glue2_msm_claude.o "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" > glue_msm_Nf${NF}_run_claude.log 2>&1 &
done
echo "launched 6 procs; waiting (F^2 ~52 min/ens, F_12 ~33 min/ens; all parallel)..." | tee -a "$BLOG"
wait

echo "================ ALL DONE ================" | tee -a "$BLOG"
date | tee -a "$BLOG"
for tag in f2 msm; do
  for nf in 2 4 6; do
    d=$(dir_of "$nf" "$tag")
    echo "  $tag Nf$nf : $(ls "$d"/F_corr.*.h5 2>/dev/null | wc -l) F_corr.h5 files" | tee -a "$BLOG"
  done
done

# NOTE (analysis): the earlier 5-config F^2 smoke (F_corr.1..5, pre-thermalization) may
# still sit in data_Nf2_..._f2/. Filter k>=100 in analysis, or rm those 5 by hand.
