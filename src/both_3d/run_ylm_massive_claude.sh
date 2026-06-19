#!/usr/bin/env bash
# Local-current Y_lm tower on the MASSIVE (Real-m / m_F) interacting ensembles (gsq=8.0, nt128 L1, Nf=2,4,6).
# Same analysis as the massless production (run_ylm_prod_claude.sh), now at REAL valence mass = SEA mass
# (the unitary point; each ensemble carries its own mRe in the dir name).  IMAG mass (m_P) is OUT OF SCOPE
# here -- both ylm programs assert(!parity); it needs the tilde-D backward leg (see jj_corr_dilute).
#
# Two programs (real mass m_F is fully in scope -> no assert trips):
#   CONN (jj_local_ylm_conn_stoch.o): per-m vector+axial tower, t0=0 wall source, spin dilution, --mass-re m.
#   DISC (jj_local_ylm_disc_stoch.o): per-m disc one-point J^a_{lm} + sigma_PS, --disc-tblock 2 (clean disc;
#     tb=2 pushes the time-dilution self-contraction bias to dt=64 only -- see disc_dilution_bias_claude.md).
#
# ORDER (user): loop ensembles, conn THEN disc for EACH so heavy masses are analyzable first.
#   Outer loop = mass, HEAVIER first: 0.200000 -> 0.100000 -> 0.050000 -> 0.010000.
#   Per mass: conn(Nf2 GPU0 ; Nf4,Nf6 GPU1) in parallel -> wait ; then disc(same layout) -> wait.
#   So each mass block finishes COMPLETELY (conn + disc, all Nf) before the next, heaviest done first.
#
# Distribution mirrors production: Nf2 -> GPU 0 ; Nf4, Nf6 -> GPU 1 (co-resident under the running CUDA MPS
# daemon).  At most 2 procs/GPU at once (12 GB TITAN V) -- conn is the heavy one; 2 conn/GPU is proven OK.
# Configs: ckpoint_lat.k, k = KMIN, KMIN+STRIDE, ... (NCONF configs); per-config atomic .h5 + resume-skip,
# so this is safe to re-launch.  Output: data_<ens>_vmRe<m>vmIm0.000000/{corr_ylm_conn_t00_nhits1_s1,
# corr_ylm_disc_tb2_nhits1}/ .  Read back: ylm_{conn,disc}_mRe<m>_nf<nf>_claude.log .
# MPS must be up: `nvidia-cuda-mps-control -d`.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

KMIN=40
STRIDE=8
NCONF=140
KMAX=$(( KMIN + STRIDE*NCONF ))      # exclusive -> 1160 ; k=40,48,...,1152 (140 configs)
NHITS=1
DISC_TB=2                            # clean disc (interval = Nt/2 = 64 ; bias only at dt=64)

MASSES="0.200000 0.100000 0.050000 0.010000"   # REAL m, HEAVIER first

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

CONN=jj_local_ylm_conn_massive.o     # distinct binary names (no ETXTBSY vs prod/disc binaries)
DISC=jj_local_ylm_disc_massive.o

# ---- build once (early-exit on failure) ----
echo "### compile -> $CONN, $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

ENS () { echo "Nf${1}_gsq8.000000at0.200000nu01.000000mRe${2}mIm0.000000nt128L1/"; }

run_conn () {   # $1=Nf  $2=GPU  $3=mass
  local NF=$1 GPU=$2 M=$3 LOG="ylm_conn_mRe${3}_nf${1}_claude.log"
  echo "### CONN Nf=$NF mRe=$M -> GPU $GPU  ($NCONF cfg, k=$KMIN..$((KMAX-STRIDE)) step $STRIDE)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$CONN" --ens-dir "$(ENS $NF $M)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --t0 0 --spin-dilution --mass-re "$M" 2>&1 | tee -a "$LOG"
  echo "### CONN Nf=$NF mRe=$M done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

run_disc () {   # $1=Nf  $2=GPU  $3=mass
  local NF=$1 GPU=$2 M=$3 LOG="ylm_disc_mRe${3}_nf${1}_claude.log"
  echo "### DISC Nf=$NF mRe=$M -> GPU $GPU  ($NCONF cfg, tblock=$DISC_TB)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$DISC" --ens-dir "$(ENS $NF $M)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" --mass-re "$M" 2>&1 | tee -a "$LOG"
  echo "### DISC Nf=$NF mRe=$M done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

for M in $MASSES; do
  echo "############################################################"
  echo "### MASS BLOCK mRe=$M  [$(date +%F_%H:%M:%S)] ###"
  echo "############################################################"

  echo "### mRe=$M PHASE conn: Nf2 GPU0 ; Nf4, Nf6 GPU1 (MPS)  [$(date +%F_%H:%M:%S)] ###"
  run_conn 2 0 "$M" &
  run_conn 4 1 "$M" &
  run_conn 6 1 "$M" &
  wait
  echo "### mRe=$M conn phase done  [$(date +%F_%H:%M:%S)] ###"

  echo "### mRe=$M PHASE disc: Nf2 GPU0 ; Nf4, Nf6 GPU1 (MPS)  [$(date +%F_%H:%M:%S)] ###"
  run_disc 2 0 "$M" &
  run_disc 4 1 "$M" &
  run_disc 6 1 "$M" &
  wait
  echo "### mRe=$M disc phase done -- ensemble mRe=$M FULLY ANALYZABLE  [$(date +%F_%H:%M:%S)] ###"
done

echo "### ALL massive Real-m ylm runs done  [$(date +%F_%H:%M:%S)] ###"
