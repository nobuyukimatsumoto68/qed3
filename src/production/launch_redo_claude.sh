#!/usr/bin/env bash
# =============================================================================
# _claude HANDOFF -- REDO production (corrected-gsq, bug-fixed) on AFFINE only.
# Builds the L=1,2 binaries from the _fermilab block driver (absolute geometry, mimics src/both_3d),
# then LAUNCHES 2-client MPS jobs (two streams per A100) via run_wrapper_redo_claude.sh. First try: 4h.
#   MASSLESS Nf{2,4,6}:  L1 gsq{0.5,1.0,1.5} (KMAX 2000)   L2 gsq{1.0,2.0,3.0} (KMAX 1200)
#   MASSIVE  Nf2      :  L1 gsq1.5 m{0.1..0.4} (KMAX 120)  L2 gsq3.0 m{0.1..0.4} (KMAX 80)
#   -> 26 streams -> 13 pair-jobs, all in /lustre2/affine/redo.
# Params: /project/qed3/qed3/src/production/params_L1L2_claude.md , params_massive_claude.md .
#
# Usage:  bash ~/launch_redo_claude.sh                    # DRY-RUN: print plan, build/submit nothing
#         bash ~/launch_redo_claude.sh --smoke            # build + submit 2 slots (1 massless + 1 massive)
#         bash ~/launch_redo_claude.sh --apply [NCHAIN]   # build + submit all 13 slots (NCHAIN default 1)
#         bash ~/launch_redo_claude.sh --rebuild --apply  # force recompile even if up-to-date
# Read back via ~/launch_redo_claude.log.
# =============================================================================
set -u
MODE=dry; NCHAIN=1; REBUILD=0
for a in "$@"; do
  case "$a" in
    --smoke)   MODE=smoke ;;
    --apply)   MODE=apply ;;
    --rebuild) REBUILD=1 ;;
    [0-9]*)    NCHAIN="$a" ;;
  esac
done

SCRIPTDIR=/project/qed3/qed3/src/production   # git-synced: driver, run/wrapper scripts, docs live here
SRC=$SCRIPTDIR/hmc_hasenbusch_block_fermilab_claude.cu
OUTDIR=/lustre2/affine/redo                   # off-git: binaries, checkpoints, logs
BINDIR=$OUTDIR
WRAP=$SCRIPTDIR/run_wrapper_redo_claude.sh
RUNSCRIPT=$SCRIPTDIR/run_redo_mps2_claude.sh
LOG=$HOME/launch_redo_claude.log
ENVSH=/home/nmatsum/env.sh; [ -f "$ENVSH" ] || ENVSH=/lustre2/affine/env.sh

# build matrix: "OUT LREF NF KMAX KRNG"
#   massless (mRe set to 0 at runtime): full trajectory targets, thin rng by 5
#   massive  (mRe set per stream)     : short targets, FULL rng checkpointing (KRNG=1)
BUILDS=(
  "hmc_fermilab_redo_massless_L1_Nf2_claude.o 1 2 2000 5"
  "hmc_fermilab_redo_massless_L1_Nf4_claude.o 1 4 2000 5"
  "hmc_fermilab_redo_massless_L1_Nf6_claude.o 1 6 2000 5"
  "hmc_fermilab_redo_massless_L2_Nf2_claude.o 2 2 1200 5"
  "hmc_fermilab_redo_massless_L2_Nf4_claude.o 2 4 1200 5"
  "hmc_fermilab_redo_massless_L2_Nf6_claude.o 2 6 1200 5"
  "hmc_fermilab_redo_massive_L1_claude.o      1 2  120 1"
  "hmc_fermilab_redo_massive_L2_claude.o      2 2   80 1"
)

{
echo "######## REDO L=1,2 build + launch ($MODE)  $(date) ########"; hostname
echo "# affine only ; 4h/job ; 2 streams/GPU (MPS) ; massless Nf{2,4,6} + massive Nf2 ; NCHAIN=$NCHAIN"
echo "# SRC=$SRC"
echo "# geometry (in _fermilab driver): /project/affine/nmatsum/qed3/geometry/  (absolute, mimics both_3d)"

# --- sanity: driver uses absolute geometry (must NOT be relative on FNAL) ---
echo; echo "== sanity: geometry paths in the _fermilab driver =="
grep -nE 'const std::string dir|#include ".*geodesic' "$SRC" | grep -vE '^\s*//' | sed 's/^/  /'

# --- build ---
echo; echo "== [1] build ${#BUILDS[@]} binaries into $BINDIR (arch sm_80) =="
if [ "$MODE" != dry ]; then
  source "$ENVSH"
  command -v nvcc >/dev/null || { echo "ERROR: nvcc not on PATH after sourcing $ENVSH -- ABORT"; exit 1; }
  : "${CUDA_PATH:=/srv/software/el8/x86_64/hpc/cuda/12.2.1}"
  NVCCFLAGS="-arch=sm_80 -g -O3 -std=c++20 -L${CUDA_PATH}/lib64 -L${CUDA_PATH}/lib64/stubs -lcudart -lcuda -lnvidia-ml -lcublas -lcufft -ldl -lgomp -Xcompiler -fopenmp"
  INCLUDES="-I/project/qed3/qed3/src/production/includes/ -I/project/qed3/qed3/qfe_mod/include/ -I/project/qed3/highfive/include/ -I/srv/software/el8/x86_64/eb/HDF5/1.14.2-GCC-12.3.0-serial/include/ -I/project/qed3/gsl/include/"
  LDFLAGS="-L/srv/software/el8/x86_64/eb/HDF5/1.14.2-GCC-12.3.0-serial/lib/ -L/project/qed3/gsl/lib/ -lhdf5 -lgsl -lgslcblas -lm"
  mkdir -p "$BINDIR"
  for spec in "${BUILDS[@]}"; do
    read -r out L NF KMAX KRNG <<< "$spec"
    dst=$BINDIR/$out
    # Incremental: (re)build only if missing, older than SRC, or --rebuild. NO rm -- nvcc overwrites in place.
    if [ "$REBUILD" -eq 0 ] && [ -x "$dst" ] && [ "$dst" -nt "$SRC" ]; then
      echo "  up-to-date: $out (LREF=$L NF=$NF KMAX=$KMAX KRNG=$KRNG)"; continue
    fi
    echo "  build $out (LREF=$L NF=$NF KMAX=$KMAX KRNG=$KRNG)  [$(date +%H:%M:%S)]"
    nvcc $NVCCFLAGS $INCLUDES -DLREF=$L -DNF=$NF -DKMAX=$KMAX -DKRNG=$KRNG "$SRC" $LDFLAGS -o "$dst" \
      2> "$BINDIR/build_${out%.o}_claude.log" \
      || { echo "  BUILD FAILED: $out (see build_${out%.o}_claude.log) -- ABORT"; exit 1; }
    [ -x "$dst" ] || { echo "  MISSING after build: $out -- ABORT"; exit 1; }
  done
  echo "== BUILD OK =="; ls -la "$BINDIR"/hmc_fermilab_redo_*.o | sed 's/^/  /'
else
  echo "  would build (incremental unless --rebuild):"
  for spec in "${BUILDS[@]}"; do read -r out L NF KMAX KRNG <<< "$spec"; echo "    $out : -DLREF=$L -DNF=$NF -DKMAX=$KMAX -DKRNG=$KRNG"; done
fi

# --- stage + permissions ---
chmod +x "$RUNSCRIPT" "$WRAP" 2>/dev/null || true

# --- launch ---
echo; echo "== [2] launch via wrapper =="
if [ "$MODE" = smoke ]; then
  echo "  SMOKE: slot 0 (massless L1 Nf2 gsq0.5+1.0), NCHAIN=1, qos=test 30min. (All massive slots 9-12 are"
  echo "         DONE/commented; pass SLOTS=... to smoke a specific active slot.) Check delta/admiss/acc, then --apply."
  QOS=test WALL=00:30:00 WSEC=1800 SLOTS="${SLOTS:-0}" bash "$WRAP" 1
elif [ "$MODE" = apply ]; then
  bash "$WRAP" "$NCHAIN"
else
  echo "  would: bash $WRAP $NCHAIN   (13 pair-slots x NCHAIN, afterany-chained per slot)"
fi

if [ "$MODE" != dry ]; then
  echo; echo "== queue =="; squeue -u "${USER:-nmatsum}" -o "%.10i %.14j %.8T %.10M %.12E" | awk 'NR==1 || $2 ~ /hmc_redo_/'
fi
echo; echo "# CHECK per-stream logs in $OUTDIR (hmc_redo_p*_L*_Nf*_gsq*_mRe*_jid*_claude.log):"
echo "#   '# gsq =', '# mass =', '# delta' small (~1.5e-5), admissibility PASS, '# alat', dH / acceptance."
echo "######## DONE ($MODE) ########"
} 2>&1 | tee "$LOG"
