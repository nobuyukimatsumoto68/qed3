#!/usr/bin/env bash
# tmp_eig_dw_irl_L2L4_claude.sh
# _claude: low complex spectrum of the BARE Wilson D_W at L2 and L4, gsq=8, thermalized -- via the IRL
# (env DW_LOWMODES): Hermitian IRL on D_W^dag D_W (CSR matvecs only, no Zolotarev -> L4 feasible) for the
# smallest-|D_W| modes, then Rayleigh-Ritz B=V^dag D_W V -> complex D_W eigenvalues near the origin.  Plotted
# with the domain wall M=1: the physics question is whether the physical branch crosses M=1 (light modes) at
# finer L, i.e. whether the overlap gets a genuine near-zero mode.  Dense geev is infeasible here
# (L2 N=10752, L4 N=41472).  GPU 1.  Reads back: eig_dw_irl_L2L4_claude.log + eig_dw_irl_L{2,4}_claude.png .
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GPU=1
NT=128
#     L | thermalized config k
LREFS=(2 4)
KS=(740 150)

# IRL params: alpha auto (power iteration), beta=2 (|D_W|<2), even deg.  Nk PER L: only ~17 modes have
# |D_W|<2 at L2 (17 converged with Nk=64), so ask for fewer than that -> all converge cleanly.
DEG=8
BETA=2.0
#        L2  L4
NKS=(16 64)
NMS=(64 192)

NVCC=nvcc
NVCCBASE="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

LOG=eig_dw_irl_L2L4_claude.log

echo "### D_W low-mode IRL  gsq=8  L=2,4  Nt=${NT}  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"

for i in "${!LREFS[@]}"; do
  L="${LREFS[$i]}"
  k="${KS[$i]}"
  NK="${NKS[$i]}"
  NM="${NMS[$i]}"
  cfg="Nf2_gsq8.000000at0.200000nu01.000000nt128L${L}/ckpoint_lat.${k}"
  N=$((2*(10*L*L+2)*NT))
  BIN="eig_arnoldi_L${L}_nt${NT}.o"
  echo "### BUILD L${L} (N=${N}, Arnoldi driver) -> $BIN ###" | tee -a "$LOG"
  $NVCC $NVCCBASE -DLREF=${L} -DNT_CLI=${NT} $INCLUDES $LDFLAGS eig_arnoldi_claude.cu -o "$BIN" 2>>"$LOG" \
    || { echo "### L${L} BUILD FAILED (see $LOG) ###" | tee -a "$LOG"; exit 1; }
  echo "### RUN L${L}  config k=${k}  m(=Nm)=${NM}  GPU$GPU  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  CONFIG_LAT="$cfg" DW_ARNOLDI=1 CUDA_VISIBLE_DEVICES=$GPU ./"$BIN" 0.0 0.0 -1 ${BETA} ${DEG} ${NK} ${NM} 2>&1 | tee -a "$LOG" \
    || { echo "### L${L} RUN FAILED ###" | tee -a "$LOG"; exit 1; }
  echo "### PLOT L${L} (bare D_W, wall M=1)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
  python3 eig_spectrum_scatter_claude.py \
    "eig_dw_arnoldi_L${L}_nt${NT}_claude.dat" \
    "eig_dw_arnoldi_L${L}_claude.png" \
    "bare D_W low modes (L${L}, gsq8, shift-invert Arnoldi)" \
    1.0 2>&1 | tee -a "$LOG"
done

echo "### D_W low-mode IRL done  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
