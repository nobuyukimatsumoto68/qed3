#!/bin/bash -l
set -u

# ============================================================================
# INTERACTING re-measurement of the shape-basis glueball correlators (F linear +
# F^2 = 0++) with the PRODUCTION 4-shape basis {triangle, rect, figure-8, three-
# triangle}, face_sign ON, ell=0..2 (n_lm=9, nops=36), TMAX_CORR=16.
#
# Full sweep = 26 ensembles (nu0=1.0, Nt=128, stride=1):
#   gsq8 x {L1,L2,L4} x Nf{2,4,6}                                     (9)
#   L1 gsq scan: Nf2 {1,2,2.4,2.5,4,12} Nf4 {1,2,2.2,2.5,4,12} Nf6 {1,2,2.4,4,12}  (17)
#
# L (N_REFINE) is a COMPILE-TIME flag -> six binaries (F/F^2 x L1/L2/L4).
# Runtime args per driver: gsq Nf nu0 kmax kmin stride.  Output h5 -> data_<cfgdir>/.
#
# ------------------------------------------------------------------------------
# PREREQUISITE (run YOURSELF -- no rm is ever placed in a script):
#   the existing glue_{msm,f2}_shapes.*.h5 in these data_ dirs are the OLD mixed
#   basis (5-shape, ell=0..3, tmax=10) and the resume-skip keys on the "complete"
#   flag, so they MUST be removed first or nothing new is measured. Remove ONLY the
#   two production prefixes (leave glue_msm_shapes_nofs.* alone):
#
#     while read cd; do
#       dd="data_$cd"
#       rm -f "$dd"/glue_msm_shapes.[0-9]*.h5 "$dd"/glue_f2_shapes.[0-9]*.h5
#     done < interacting_ens_list_claude.txt
#
#   The python3 preflight below ABORTS if any stale old-basis h5 (missing the
#   n_lm attribute) is still present, so a forgotten rm cannot silently no-op.
# ------------------------------------------------------------------------------
# Reference: free-limit anchors F ell=1 -> sqrt(2)=1.41421 (L=1 det 1.33242);
#            F^2 0++ -> two-photon 2sqrt(2)=2.828 (see run_free_validation_claude.sh).
# ============================================================================

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
# CPU-only job (gradient flow + host Wilson-loop holonomies; no GPU kernels). Run 8 ensembles
# concurrently, ONE single-threaded process each => 8 cores total.
export OMP_NUM_THREADS=1
NWORK=8

NU0=1.0
NU0STR=1.000000    # std::to_string(1.0) form used in the dir names
KMAX=100000
KMIN=1
STRIDE=1

NVCC=/usr/local/cuda-12.6/bin/nvcc
NVCCFLAGS="-w -arch=sm_70 -O2 -lcusolver -std=c++17 -Xcompiler -fopenmp"
H5I="-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
INCLUDES="-Iincludes -I../../qfe_mod/include -I/usr/local/cuda-12.6/include/ $H5I"
LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/12.3/targets/x86_64-linux/lib/ -L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5"

LOG=interacting_shapes_claude.log
echo "================ INTERACTING shapes sweep START $(date) ================" | tee "$LOG"

# ---- ensemble table: "Nf gsq L" (nu0=1.0), ordered LARGEST-nckpt-FIRST so the 8-worker pool
#      grabs the long jobs early (better load balance; no single giant job finishes alone last).
#      Trailing comment = approx nckpt at survey time.
# (Also written to interacting_ens_list_claude.txt in cfgdir form for the rm above.)
ENS=(
  "4 2.000000 1"  # 5899
  "6 2.000000 1"  # 5538
  "2 2.000000 1"  # 4440
  "4 8.000000 1"  # 3557
  "2 8.000000 1"  # 3426
  "6 8.000000 1"  # 3163
  "2 2.500000 1"  # 2723
  "4 2.200000 1"  # 2092
  "2 8.000000 2"  # 1601
  "2 2.400000 1"  # 1418
  "6 2.400000 1"  # 1264
  "4 2.500000 1"  # 1099
  "4 8.000000 2"  # 518
  "2 1.000000 1"  # 319
  "2 12.000000 1" # 319
  "2 4.000000 1"  # 319
  "4 1.000000 1"  # 319
  "4 12.000000 1" # 319
  "6 1.000000 1"  # 319
  "6 12.000000 1" # 319
  "6 8.000000 2"  # 289
  "4 8.000000 4"  # 190
  "4 4.000000 1"  # 184
  "2 8.000000 4"  # 168
  "6 4.000000 1"  # 153
  "6 8.000000 4"  # 68
)

# cfgdir list for the manual rm prerequisite
: > interacting_ens_list_claude.txt
for e in "${ENS[@]}"; do
  set -- $e
  echo "Nf${1}_gsq${2}at0.200000nu0${NU0STR}nt128L${3}" >> interacting_ens_list_claude.txt
done

# ---- build six binaries -----------------------------------------------------
echo "---- build $(date) ----" | tee -a "$LOG"
for L in 1 2 4; do
  if [ "$L" = "1" ]; then DEF=""; else DEF="-DN_REFINE_CLI=$L"; fi
  "$NVCC" glue2_msm_shapes_claude.cu $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue2_msm_shapes_L${L}_claude.o 2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F  L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
  "$NVCC" glue_f2_shapes_claude.cu   $DEF $NVCCFLAGS $INCLUDES $LDFLAGS -o glue_f2_shapes_L${L}_claude.o   2>&1 | tee -a "$LOG"
  [ "${PIPESTATUS[0]}" -ne 0 ] && { echo "F^2 L$L BUILD FAILED" | tee -a "$LOG"; exit 1; }
done

# ---- preflight: abort if any stale old-basis h5 (missing n_lm attr) remains --
echo "---- preflight (stale old-basis h5 check) $(date) ----" | tee -a "$LOG"
python3 - "$NU0STR" <<'PY' 2>&1 | tee -a "$LOG"
import sys, glob, os
import h5py
nu0 = sys.argv[1]
# rebuild cfgdirs from the same table the shell uses
table = [
  (2,"8.000000",1),(4,"8.000000",1),(6,"8.000000",1),
  (2,"8.000000",2),(4,"8.000000",2),(6,"8.000000",2),
  (2,"8.000000",4),(4,"8.000000",4),(6,"8.000000",4),
  (2,"1.000000",1),(2,"2.000000",1),(2,"2.400000",1),(2,"2.500000",1),(2,"4.000000",1),(2,"12.000000",1),
  (4,"1.000000",1),(4,"2.000000",1),(4,"2.200000",1),(4,"2.500000",1),(4,"4.000000",1),(4,"12.000000",1),
  (6,"1.000000",1),(6,"2.000000",1),(6,"2.400000",1),(6,"4.000000",1),(6,"12.000000",1),
]
stale = []
for nf,g,L in table:
    dd = f"data_Nf{nf}_gsq{g}at0.200000nu0{nu0}nt128L{L}"
    for pref in ("glue_msm_shapes","glue_f2_shapes"):
        fs = sorted(glob.glob(f"{dd}/{pref}.[0-9]*.h5"))
        if not fs:
            continue
        try:
            with h5py.File(fs[0],"r") as h:
                has_nlm = ("n_lm" in h) or ("n_lm" in h.attrs)
        except Exception as ex:
            has_nlm = False
        if not has_nlm:
            stale.append(f"{dd}/{pref}.*.h5  (sample {os.path.basename(fs[0])} lacks n_lm; {len(fs)} files)")
if stale:
    print("PREFLIGHT ABORT: stale OLD-basis h5 still present. Remove these first (run yourself):")
    for s in stale:
        print("   ", s)
    sys.exit(3)
print("preflight OK: no stale old-basis h5 detected.")
PY
[ "${PIPESTATUS[0]}" -ne 0 ] && { echo "PREFLIGHT FAILED -- remove stale h5 and rerun." | tee -a "$LOG"; exit 1; }

# ---- measure: 8-worker pool, one single-threaded process per ensemble --------
# Each worker owns an ensemble: runs F then F^2 for it, writes its own per-ensemble
# log, then xargs hands it the next ensemble. Master LOG records start/done stamps.
run_one() {
  local NF=$1
  local GSQ=$2
  local L=$3
  local tag="Nf${NF}_gsq${GSQ}_L${L}"
  local elog="interacting_shapes_${tag}_claude.log"
  echo "[start $(date '+%F %T')] $tag" >> "$LOG"
  {
    echo "==== $tag  F (linear)  $(date) ===="
    ./glue2_msm_shapes_L${L}_claude.o "$GSQ" "$NF" "$NU0" "$KMAX" "$KMIN" "$STRIDE"
    echo "==== $tag  F^2 (squared)  $(date) ===="
    ./glue_f2_shapes_L${L}_claude.o   "$GSQ" "$NF" "$NU0" "$KMAX" "$KMIN" "$STRIDE"
  } > "$elog" 2>&1
  echo "[done  $(date '+%F %T')] $tag  (log: $elog)" >> "$LOG"
}
export -f run_one
export NU0 KMAX KMIN STRIDE LOG

echo "---- launch $NWORK-worker pool $(date) ----" | tee -a "$LOG"
printf '%s\n' "${ENS[@]}" | xargs -P "$NWORK" -I{} bash -c 'run_one $1' _ {}

echo "================ INTERACTING shapes sweep DONE $(date) ================" | tee -a "$LOG"
