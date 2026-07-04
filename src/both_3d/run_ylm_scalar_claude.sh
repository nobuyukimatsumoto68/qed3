#!/usr/bin/env bash
# PRODUCTION scalar-density (sigma_PS, sigma_FS) Y_lm tower, alongside / appended to the local-current tower.
# Mirrors run_ylm_prod_claude.sh but uses the SCALAR drivers:
#   CONN (jj_local_ylm_scalar_conn_stoch.o): vector+axial tower AND the a=0 scalar V++ (h0/scalar/...).
#   DISC (jj_local_ylm_scalar_disc_stoch.o): disc J^a + sigma_PS condensate AND the a=0 complex loop
#                                            J^0_{lm} (h0/disc/ylm/s0/...).
# Physical scalar correlator (analysis) = (N_f/2)(-C_c + C_d);  C_c = -(Vpp+Vmm) from h0/scalar,
#   C_d = two_point(J^0);  Re(J^0)->sigma_PS, Im(J^0)->sigma_FS.  Massless / m_F only (m_P out of scope).
#
# MODE (positional $1, default "together"):
#   together   -- compute vector+axial+scalar FRESH in one file (for ensembles WITHOUT the jj tower yet).
#                 Output = the SAME dirs the original conn/disc drivers use (superset), with "complete".
#   scalaronly -- add ONLY the scalar keys and APPEND into the EXISTING conn/disc per-hit .h5 (for ensembles
#                 that ALREADY have the jj tower).  Must match the original run's t0/spin/nhits/disc-tblock so
#                 the output paths + RNG seeds line up (the wall/volume source is regenerated identically).
#
# Ensemble knobs (env-overridable): GSQ, NREF (1 or 2; sets -DN_REFINE_CLI + the L tag), KMIN/STRIDE/NCONF.
# Distribution: Nf2 -> GPU 0 ; Nf4, Nf6 -> GPU 1 (<=2 procs/GPU under the running CUDA MPS daemon).
# Read back: ylm_scalar_{conn,disc}_nf{2,4,6}_claude.log .   (MPS up: `nvidia-cuda-mps-control -d`.)
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

MODE="${1:-together}"
case "$MODE" in
  together)   SCALAR_FLAG="" ;;
  scalaronly) SCALAR_FLAG="--IsScalarOnly" ;;
  *) echo "usage: $0 [together|scalaronly]"; exit 1 ;;
esac

GSQ="${GSQ:-8.000000}"
NREF="${NREF:-1}"
LTAG="L${NREF}"
KMIN="${KMIN:-40}"; STRIDE="${STRIDE:-8}"; NCONF="${NCONF:-140}"
KMAX=$(( KMIN + STRIDE*NCONF ))     # exclusive
NHITS="${NHITS:-1}"
# DISC_TB MUST match the existing disc dir (corr_ylm_disc_tb<TB>_nhits<H>) for scalaronly APPEND to land in
# the right file.  The existing ylm production disc data is tb2 -> default 2.  (together mode: any tb is fine.)
DISC_TB="${DISC_TB:-2}"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -DN_REFINE_CLI=${NREF} -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

CONN=jj_local_ylm_scalar_conn_stoch.o
DISC=jj_local_ylm_scalar_disc_stoch.o

echo "### MODE=$MODE  GSQ=$GSQ  $LTAG  k=$KMIN..$((KMAX-STRIDE)) step $STRIDE ($NCONF cfg)  [$(date +%F_%H:%M:%S)] ###"
echo "### compile -> $CONN, $DISC (NREF=$NREF) ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_scalar_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_scalar_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

ENS () { echo "Nf${1}_gsq${GSQ}at0.200000nu01.000000nt128${LTAG}/"; }

run_conn () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2 LOG="ylm_scalar_conn_nf${1}_claude.log"
  echo "### CONN Nf=$NF -> GPU $GPU ($MODE)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$CONN" --ens-dir "$(ENS $NF)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --t0 0 --spin-dilution $SCALAR_FLAG 2>&1 | tee -a "$LOG"
  echo "### CONN Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

run_disc () {   # $1=Nf  $2=GPU
  local NF=$1 GPU=$2 LOG="ylm_scalar_disc_nf${1}_claude.log"
  echo "### DISC Nf=$NF -> GPU $GPU ($MODE)  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
  CUDA_VISIBLE_DEVICES=$GPU ./"$DISC" --ens-dir "$(ENS $NF)" --kmin "$KMIN" --stride "$STRIDE" --kmax "$KMAX" \
    --nhits "$NHITS" --disc-tblock "$DISC_TB" $SCALAR_FLAG 2>&1 | tee -a "$LOG"
  echo "### DISC Nf=$NF done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
}

echo "### PHASE 1: CONN  [$(date +%F_%H:%M:%S)] ###"
run_conn 2 0 &
run_conn 4 1 &
run_conn 6 1 &
wait
echo "### PHASE 1 (conn) done  [$(date +%F_%H:%M:%S)] ###"

echo "### PHASE 2: DISC  [$(date +%F_%H:%M:%S)] ###"
run_disc 2 0 &
run_disc 4 1 &
run_disc 6 1 &
wait
echo "### ALL scalar ylm ($MODE) done  [$(date +%F_%H:%M:%S)] ###"
