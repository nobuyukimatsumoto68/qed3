#!/bin/bash -l
set -u

# ============================================================================
# FREE-ENSEMBLE shapes measurement for BOTH F (linear F_12) and F^2 (0++), with the
# 4-shape untwisted basis {triangle, rect, figure-8, three-triangle}, face_sign ON,
# tmax=16, on the free (non-compact Maxwell) ensemble gsq8...nt128L1 (Nf=0).
#
# stride=10 -> all 10000 configs (k=10,20,...,100000). Output into data_gsq8...nt128L1/:
#   glue_msm_shapes.<k>.h5   (F,  linear)
#   glue_f2_shapes.<k>.h5    (F^2, squared)
#
# PREREQUISITE (run yourself -- no rm in scripts): the OLD free h5 are a DIFFERENT basis
# (3-shape, and/or tmax=10) and the resume-skip keys on the complete flag, so remove them first:
#     rm data_gsq8.000000at0.200000nt128L1/glue_msm_shapes.*.h5
#     rm data_gsq8.000000at0.200000nt128L1/glue_msm_shapes_nofs.*.h5
#   (glue_f2_shapes.* on the free dir do not exist yet, nothing to remove there.)
#
# Reference: F l=1 -> sqrt(2)=1.41421 (paper Delta_0; L=1 lattice ~1.33242).
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=1

# ---- knobs ----
GSQ=8.000000
NF=0
NU0=1.0
KMIN=10
STRIDE=10       # smallest = config spacing -> ALL 10000 configs
KMAXRUN=100000
FREEDIR="gsq8.000000at0.200000nt128L1"
DATADIR="data_${FREEDIR}"
# analysis
AT=0.2
TCUT=14
BINSIZE=0
RTOL=1e-8
TREBASE=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=free_validation_claude.log
echo "================ FREE F + F^2 (4-shape) START $(date) ================" | tee "$LOG"

# ---- build (distinct binary names; do not touch the running sweep's binaries) ----
g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o glue_gevp_analysis_claude.o 2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "analysis BUILD FAILED" | tee -a "$LOG"; exit 1; }
"$NVCC" glue2_msm_shapes_claude.cu $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_freeON_claude.o 2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F  driver BUILD FAILED" | tee -a "$LOG"; exit 1; }
"$NVCC" glue_f2_shapes_claude.cu   $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_freeON_claude.o   2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F^2 driver BUILD FAILED" | tee -a "$LOG"; exit 1; }

# ---- measure both on the free ensemble (Nf=0) ----
echo "================ MEASURE F (linear) $(date) ================" | tee -a "$LOG"
./glue2_msm_shapes_freeON_claude.o "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" 2>&1 | tail -2 | tee -a "$LOG"
echo "================ MEASURE F^2 (squared) $(date) ================" | tee -a "$LOG"
./glue_f2_shapes_freeON_claude.o   "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" 2>&1 | tail -2 | tee -a "$LOG"

# ---- analyze ----
# arg order: dir NF NCH drop_l0 at outfile nops2 tcut binsize kmin rtol prefix trebase orbits vacsub keepl0 evec_t lsel use_rebase
echo "================ ANALYZE $(date) ================" | tee -a "$LOG"
# NO vacuum subtraction anywhere (vacsub=0): a no-op for F (l>=1, <O>=0), and for F^2 the GEVP
# separates the constant/vacuum mode on its own (identical signal + errors, simpler). Args tail:
#   ... prefix trebase orbits vacsub keepl0 evec_t lsel use_rebase omit_top
# F: l=1 sector, rebase trebase=1 -> 3 physical states (the m-triplet); no constant mode -> omit_top=0
./glue_gevp_analysis_claude.o "$DATADIR" 1 1 0 "$AT" gevp_msm_free_l1_n3_claude.dat 3 "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" glue_msm_shapes "$TREBASE" "" 0 1 0 "1" 1 0 2>&1 | tee -a "$LOG"
printf "F  l=1 3-state : " | tee -a "$LOG"; awk 'NR>1{printf "%.3f(%.3f) ",$2,$3}' gevp_msm_free_l1_n3_claude.dat | tee -a "$LOG"; echo "" | tee -a "$LOG"
# F^2: 0++ scalar = l=0 sector; rebase to 2 states = constant/vacuum mode (ground, Delta~0) + the
# lightest PHYSICAL 0++ (state 0 = two-photon 2sqrt2). nops2=2, omit_top=0.
./glue_gevp_analysis_claude.o "$DATADIR" 1 1 0 "$AT" gevp_f2_free_l0_claude.dat  2 "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" glue_f2_shapes  "$TREBASE" "" 0 1 0 "0" 1 0 2>&1 | tee -a "$LOG"
printf "F^2 0++ (2sqrt2): " | tee -a "$LOG"; awk 'NR>1{printf "%.3f(%.3f) ",$4,$5}' gevp_f2_free_l0_claude.dat  | tee -a "$LOG"; echo "" | tee -a "$LOG"
echo "REFERENCE: F l=1 -> sqrt(2)=1.41421 (L=1 det 1.33242)" | tee -a "$LOG"
echo "================ DONE $(date) ================" | tee -a "$LOG"
