#!/bin/bash -l
set -u

# ============================================================================
# REPRESENTATIVE-FIRST shape-basis glueball run.
#
# Focus ensemble (foreground): linear F_12 (msm), Nf=2, L=1, gsq=8.
#   stride=1 (ALL configs), kmax as large as possible (=KMAXRUN), kmin=KMIN.
#   -> measure -> GEVP -> gevp_msm_Nf2_gsq8.000000_L1_claude.dat  (eyeball this first).
# Then the remaining ensembles (both ops, all massless L=1,2,4, all gsq) are launched in the
# BACKGROUND via run_glue_shapes_all_claude.sh (resume-safe; skips the representative already done).
#
# PREREQUISITE (run yourself -- no rm in scripts, per project rule): the OLD shapes h5 are
# INCOMPATIBLE now (linear driver gained ell=0 -> n_lm 15->16; both drivers window dt=0..TMAX).
# Remove them first, e.g.:
#     rm data_Nf*_gsq*at0.200000nu0*nt128L*/glue_msm_shapes.*.h5
#     rm data_Nf*_gsq*at0.200000nu0*nt128L*/glue_f2_shapes.*.h5
# (The mesonic correlators and any F_corr.* density files are untouched.)
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=1

# ---- knobs ----
STRIDE=1
KMIN=20
KMAXRUN=1000000        # "as large as possible": >> max config index, so ALL configs are used
NU0=1.0
# analysis
AT=0.2
NOPS2=8
TCUT=5
BINSIZE=0
RTOL=1e-8

# representative ensemble
NF=2
GSQ=8.000000
REP="Nf${NF}_gsq${GSQ}at0.200000nu01.000000nt128L1"   # nu0 = 1.0 -> "nu01.000000" (6 dp, matches driver)

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=glue_shapes_repr_claude.log
echo "================ REPR START $(date) ================" | tee "$LOG"

# ---- build: analysis (g++) + L=1 linear driver (nvcc) ----
g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o glue_gevp_analysis_claude.o 2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "analysis BUILD FAILED" | tee -a "$LOG"; exit 1; }
"$NVCC" glue2_msm_shapes_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_claude.o 2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "L=1 linear BUILD FAILED" | tee -a "$LOG"; exit 1; }

# ---- 1) representative: measure (foreground) then analyze ----
echo "================ REPR MEASURE msm $REP $(date) ================" | tee -a "$LOG"
./glue2_msm_shapes_claude.o "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" 2>&1 | tee -a "$LOG"
echo "================ REPR ANALYZE $(date) ================" | tee -a "$LOG"
./glue_gevp_analysis_claude.o "data_$REP" 1 1 0 "$AT" "gevp_msm_Nf${NF}_gsq${GSQ}_L1_claude.dat" "$NOPS2" "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" glue_msm_shapes 2>&1 | tee -a "$LOG"
echo "REPR .dat:" | tee -a "$LOG"
head -5 "gevp_msm_Nf${NF}_gsq${GSQ}_L1_claude.dat" | tee -a "$LOG"

# ---- 2) background: the rest (both ops, all L/gsq), stride=1, resume-safe ----
echo "================ LAUNCH BACKGROUND (rest) $(date) ================" | tee -a "$LOG"
STRIDE="$STRIDE" KMIN="$KMIN" nohup bash run_glue_shapes_all_claude.sh > glue_shapes_all_bg_claude.log 2>&1 &
echo "background PID $! -> glue_shapes_all_bg_claude.log" | tee -a "$LOG"
echo "================ REPR DONE $(date) ================" | tee -a "$LOG"
