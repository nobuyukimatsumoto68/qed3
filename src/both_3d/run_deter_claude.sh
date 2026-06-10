#!/usr/bin/env bash
# run_deter_claude.sh -- orchestrate the EXACT current-correlator validation (free-field conformal).
# Three stages (jj_exact_freefield_impl_plan_claude.md), one binary per L (-DN_REFINE_CLI):
#   1a  jj_propagator_deter   -> D_m^{-1}  (LU, mass-DEPENDENT)   data_<ESNID>/prop_deter_L<L>/
#   1b  jj_kbuild_exact       -> K^P(t)    (overlap, mass-FREE)   data_<TAG>_Kdense_L<L>/   (built once, reused)
#   2   jj_contract_deter     -> corr      (pure cuBLAS)          data_<ESNID>/corr_deter_exact_L<L>/
# K (1b) is mass-independent => built once per (config,L); propagator (1a) + contract (2) run per mass.
#
# Usage:
#   bash run_deter_claude.sh [--L "1"] [--masses "0:0"] [--ens-dir DIR|--free] [--gpu 0]
#                            [--n-t0 2] [--ninter 10] [--no-build]
#   --masses : space-separated re:im pairs, e.g. "0:0 0.1:0 0:0.1"
#   --free   : free field (default; omit --ens-dir).  --ens-dir DIR : interacting ensemble.
set -u

LS="1"; MASSES="0:0"; ENS_DIR=""; GPU=0; NT0=2; NINTER=10; DO_BUILD=1; MODE="exact"
while [[ $# -gt 0 ]]; do case "$1" in
  --L)        LS="$2"; shift 2;;
  --masses)   MASSES="$2"; shift 2;;
  --ens-dir)  ENS_DIR="$2"; shift 2;;
  --free)     ENS_DIR=""; shift;;
  --gpu)      GPU="$2"; shift 2;;
  --n-t0)     NT0="$2"; shift 2;;
  --ninter)   NINTER="$2"; shift 2;;
  --local)    MODE="local"; shift;;   # ultralocal current (no exact K); propagator + jj_local_deter
  --exact)    MODE="exact"; shift;;   # exact overlap current (default): kbuild + propagator + contract
  --no-build) DO_BUILD=0; shift;;
  -h|--help)  sed -n '2,16p' "$0"; exit 0;;
  *) echo "unknown arg: $1"; exit 1;;
esac; done

cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
module load cuda/12.8 2>/dev/null; module load gcc/13.2.0 2>/dev/null
export CUDA_VISIBLE_DEVICES="$GPU"

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

ens_args(){ [[ -z "$ENS_DIR" ]] && echo "" || echo "--ens-dir $ENS_DIR --ninter $NINTER"; }

for L in $LS; do
  KB=jj_kbuild_exact_L${L}.o
  PR=jj_propagator_deter_L${L}.o
  CT=jj_contract_deter_L${L}.o
  LC=jj_local_deter_L${L}.o
  echo; echo "############## L=${L}  (mode=${MODE}) ##############"
  # which sources to compile for this mode
  if [[ "$MODE" == "local" ]]; then SRCS="jj_propagator_deter_claude.cu:$PR jj_local_deter_claude.cu:$LC";
  else SRCS="jj_kbuild_exact_claude.cu:$KB jj_propagator_deter_claude.cu:$PR jj_contract_deter_claude.cu:$CT"; fi
  if [[ "$DO_BUILD" == 1 ]]; then
    echo "### compile (N_REFINE_CLI=${L}) ###"
    for pair in $SRCS; do
      src="${pair%%:*}"; bin="${pair##*:}"
      echo "+ nvcc ... -DN_REFINE_CLI=${L} ${src} -o ${bin}"
      $NVCC $NVCCFLAGS $INCLUDES $LDFLAGS -DN_REFINE_CLI=${L} "$src" -o "$bin" || { echo "BUILD FAILED: $src"; exit 1; }
    done
  fi

  if [[ "$MODE" == "exact" ]]; then
    echo "### [1b] exact K-build (mass-independent, once) ###"
    ./"$KB" $(ens_args) || { echo "kbuild FAILED"; exit 1; }
  fi

  for m in $MASSES; do
    MRE="${m%%:*}"; MIM="${m##*:}"
    echo "### mass (re=${MRE}, im=${MIM}) ###"
    echo "  [1a] propagator"
    ./"$PR" --mass-re "$MRE" --mass-im "$MIM" $(ens_args) || { echo "propagator FAILED"; exit 1; }
    if [[ "$MODE" == "local" ]]; then
      echo "  [2L] local-current contract"
      ./"$LC" --mass-re "$MRE" --mass-im "$MIM" --n-t0 "$NT0" $(ens_args) || { echo "local FAILED"; exit 1; }
    else
      echo "  [2] exact-current contract"
      ./"$CT" --mass-re "$MRE" --mass-im "$MIM" --n-t0 "$NT0" $(ens_args) || { echo "contract FAILED"; exit 1; }
    fi
  done
done
echo; echo "### done ###"
