#!/usr/bin/env bash
# _claude: PER-FRAME FORCE NORMS (L2/Linf) for the L4 Hasenbusch ladder -- base {0,0.4,1.0} MG20 vs safe
# {0,0.5,1.0} MG100 -- on an existing L4 config, run at Nf=6 AND Nf=2. Prints the gauge force (Nf-indep) +
# per-frame force norms (sum over Nf/2 pf stacks; ~3x at Nf6) + one MD integration -> dH / N_CG / Osborn Cost.
# Directly shows the Nf=6 concern (3 fermion forces add). Driver test_hasenbusch_forcenorm_l4_claude.cu. No rm.
#
# Args to the binary: [gsq] [Nf] [N_CFG] [N_DRAW]. NFLIST env sets the Nf sweep (default "6 2"). GSQ env.
set -u

cd /mnt/barracuda22/qed3/qed3/src/tuning

LOG=test_hasenbusch_forcenorm_l4_claude.log
GPU=${GPU:-0}
export CUDA_VISIBLE_DEVICES=$GPU
export OMP_NUM_THREADS=4
GSQ=${GSQ:-6.0}      # stiff coupling; override e.g. GSQ=2.0
# gauge config is Nf-INDEPENDENT -> load the existing Nf2 gsq<GSQ> L4 config for BOTH Nf runs
GF=$(printf "%.6f" "$GSQ")
export CFGDIR="Nf2_gsq${GF}at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000/"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_hasenbusch_forcenorm_l4_claude.cu
BIN=test_hb_forcenorm_l4.x

{
  echo "######## L4 per-frame force norms (OLD vs RETUNED ladder)  $(date)  GPU=${GPU}  gsq=${GSQ} ########"
  echo "==================== BUILD (LREF=4) ===================="
  if ! $NVCC $NVCCFLAGS $INCLUDES -DLREF=4 $SRC $LDFLAGS -o "$BIN" ; then
    echo "!!! BUILD FAILED"; exit 1
  fi
  for NF in ${NFLIST:-6 2}; do
    echo "-------------------- RUN (gsq=${GSQ} Nf=${NF} N_CFG=1 N_DRAW=1) --------------------"
    ./"$BIN" "$GSQ" "$NF" 1 1
  done
  echo "######## done ########"
} 2>&1 | tee "$LOG"
