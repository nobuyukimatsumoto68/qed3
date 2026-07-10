#!/bin/bash -l
set -u

# ============================================================================
# PHASE 3 (operator CONSOLIDATION): re-measure the shape-basis glueball correlators with the
# orbits SUMMED (equal weight per orbit) into ONE operator per shape -> 5-shape basis at every L.
# Cuts F_corr_blk from n_lm*n_orbits^2 to n_lm*25 (L2 7.8x, L4 70x smaller).
#
# ONLY L2 and L4 need re-measuring: at L1 each shape is already a single orbit, so consolidation is
# a no-op and the existing L1 phase-2 h5 (n_orbits=5) are unchanged. So the sweep = 6 ensembles:
#   gsq8 x {L2,L4} x Nf{2,4,6}   (nu0=1.0, Nt=128, stride=1)
#
# L (N_REFINE) compile-time flag -> 4 binaries (F/F^2 x L2/L4). Args: gsq Nf nu0 kmax kmin stride.
#
# ------------------------------------------------------------------------------
# PREREQUISITE (run YOURSELF -- no rm in a script): remove the phase-2 (per-orbit, n_orbits=14/42)
# h5 in the 6 L2/L4 data_ dirs first (the resume-skip keys on "complete"):
#     while read cd; do
#       dd="data_$cd"
#       rm -f "$dd"/glue_msm_shapes.[0-9]*.h5 "$dd"/glue_f2_shapes.[0-9]*.h5
#     done < consolidate_ens_list_claude.txt
#   The python3 preflight below ABORTS if any stale h5 with n_orbits != 5 remains.
#   (Do NOT touch L1 -- those are already consolidated.)
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=1
NWORK=6
NU0=1.0
NU0STR=1.000000
KMAX=100000
KMIN=1
STRIDE=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"
LOG=consolidate_remeasure_claude.log
echo "================ CONSOLIDATION re-measure (L2/L4) START $(date) ================" | tee "$LOG"

# ensembles: "Nf gsq L" (largest-nckpt first)
ENS=( "2 8.000000 2" "4 8.000000 2" "6 8.000000 2" "4 8.000000 4" "2 8.000000 4" "6 8.000000 4" )
: > consolidate_ens_list_claude.txt
for e in "${ENS[@]}"; do
  set -- $e
  echo "Nf${1}_gsq${2}at0.200000nu0${NU0STR}nt128L${3}" >> consolidate_ens_list_claude.txt
done

# ---- build L2, L4 binaries ----
echo "---- build $(date) ----" | tee -a "$LOG"
for L in 2 4; do
  "$NVCC" glue2_msm_shapes_claude.cu -DN_REFINE_CLI=$L $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_L${L}_claude.o 2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F  L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
  "$NVCC" glue_f2_shapes_claude.cu   -DN_REFINE_CLI=$L $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_L${L}_claude.o   2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F^2 L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
done

# ---- preflight: abort if any stale non-consolidated h5 (n_orbits != 5) remains ----
echo "---- preflight $(date) ----" | tee -a "$LOG"
python3 - "$NU0STR" <<'PY' 2>&1 | tee -a "$LOG"
import sys, glob, os
import h5py
nu0 = sys.argv[1]
table = [(2,"8.000000",2),(4,"8.000000",2),(6,"8.000000",2),(2,"8.000000",4),(4,"8.000000",4),(6,"8.000000",4)]
stale = []
for nf,g,L in table:
    dd = f"data_Nf{nf}_gsq{g}at0.200000nu0{nu0}nt128L{L}"
    for pref in ("glue_msm_shapes","glue_f2_shapes"):
        fs = sorted(glob.glob(f"{dd}/{pref}.[0-9]*.h5"))
        if not fs: continue
        try:
            with h5py.File(fs[0],"r") as h:
                nob = int(h["n_orbits"][0]) if "n_orbits" in h else -1
        except Exception:
            nob = -1
        if nob != 5:
            stale.append(f"{dd}/{pref}.*.h5  (sample n_orbits={nob} != 5; {len(fs)} files)")
if stale:
    print("PREFLIGHT ABORT: stale non-consolidated h5 present. Remove first:")
    for s in stale: print("   ", s)
    sys.exit(3)
print("preflight OK.")
PY
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "PREFLIGHT FAILED." | tee -a "$LOG"; exit 1; }

# ---- measure (6-worker pool) ----
run_one() {
  local NF=$1 GSQ=$2 L=$3
  local tag="Nf${NF}_gsq${GSQ}_L${L}"
  local elog="consolidate_${tag}_claude.log"
  echo "[start $(date '+%F %T')] $tag" >> "$LOG"
  {
    ./glue2_msm_shapes_L${L}_claude.o "$GSQ" "$NF" "$NU0" "$KMAX" "$KMIN" "$STRIDE"
    ./glue_f2_shapes_L${L}_claude.o   "$GSQ" "$NF" "$NU0" "$KMAX" "$KMIN" "$STRIDE"
  } > "$elog" 2>&1
  echo "[done  $(date '+%F %T')] $tag" >> "$LOG"
}
export -f run_one
export NU0 KMAX KMIN STRIDE LOG
printf '%s\n' "${ENS[@]}" | xargs -P "$NWORK" -I{} bash -c 'run_one $1' _ {}
echo "================ CONSOLIDATION re-measure DONE $(date) ================" | tee -a "$LOG"
