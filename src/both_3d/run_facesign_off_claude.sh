#!/bin/bash -l
set -u

# ============================================================================
# face_sign ISOLATION test (reproduce glue_ylms3's operator, one knob at a time).
# Measures the linear F_12 triangle operator with face_sign OFF (raw plaquette
# orientation = glue_ylms3's plaquette_angle_avg_Ylm_real), SINGLE flow=1.0, on
# the representative ensemble Nf2 L1 gsq8, then analyzes triangle-only (orbits=0).
#
# Compare the resulting ground to:
#   - face_sign ON  triangle (already have): ~1.25
#   - glue_ylms3 (density, no face_sign, 3-flow): ~1.4-1.6
# If face_sign OFF jumps to ~1.5, the orientation is the cause.
#
# Output prefix glue_msm_shapes_nofs.*.h5 (distinct -> does NOT touch production
# glue_msm_shapes.*.h5). No rm needed. NPARALLEL/OMP single.
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=1

# knobs (match the representative)
NF=2
GSQ=8.000000
NU0=1.0
STRIDE=1
KMIN=20
KMAXRUN=1000000
REP="Nf${NF}_gsq${GSQ}at0.200000nu01.000000nt128L1"
# analysis
AT=0.2
NOPS2=8
TCUT=6
BINSIZE=0
RTOL=1e-8
TREBASE=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
H5L="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=facesign_off_claude.log
echo "================ face_sign OFF test START $(date) ================" | tee "$LOG"

# build: analysis (g++) + NO_FACE_SIGN driver (nvcc)
g++ -O2 -std=c++17 -x c++ -I../../qfe_mod/include $H5I glue_gevp_analysis_claude.cu $H5L -o glue_gevp_analysis_claude.o 2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "analysis BUILD FAILED" | tee -a "$LOG"; exit 1; }
"$NVCC" glue2_msm_shapes_claude.cu $NVCCFLAGS -DNO_FACE_SIGN $INCLUDES $LDFLAGS -o glue2_msm_shapes_nofs_claude.o 2>&1 | tee -a "$LOG"
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "NO_FACE_SIGN driver BUILD FAILED" | tee -a "$LOG"; exit 1; }

# measure (single flow, face_sign OFF) -> glue_msm_shapes_nofs.*.h5
echo "================ MEASURE nofs $REP $(date) ================" | tee -a "$LOG"
./glue2_msm_shapes_nofs_claude.o "$GSQ" "$NF" "$NU0" "$KMAXRUN" "$KMIN" "$STRIDE" 2>&1 | tee -a "$LOG"

# analyze triangle-only (orbits=0), l=0..3, vacsub on
echo "================ ANALYZE nofs triangle-only $(date) ================" | tee -a "$LOG"
./glue_gevp_analysis_claude.o "data_$REP" 1 1 0 "$AT" "gevp_msm_nofs_tri_Nf2_gsq8_L1_claude.dat" "$NOPS2" "$TCUT" "$BINSIZE" "$KMIN" "$RTOL" glue_msm_shapes_nofs "$TREBASE" 0 1 1 2>&1 | tee -a "$LOG"
echo "nofs triangle ground:" | tee -a "$LOG"
awk 'NR>1{printf "%.3f(%.3f) ",$2,$3}' gevp_msm_nofs_tri_Nf2_gsq8_L1_claude.dat | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "================ DONE $(date) ================" | tee -a "$LOG"
