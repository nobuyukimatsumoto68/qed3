#!/usr/bin/env bash
# _claude: mixed-prec MULTISHIFT validation on 3 REAL thermalized configs per L (L1/L2/L4), Nf6 gsq8
# massless, Nt=128 (the configs' real temporal extent -- 8x bigger than the nt16 synthetic test, and
# with the real near-zero Wilson mode). Purpose: check whether the "shelve mixed-prec" conclusion
# (co-linearity breakdown + ~1.1x latency-bound) holds on real gauge fields, esp. L4. Reads back via
# the .log. No rm.
set -u

cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=tmp_multishift_mixed_realcfg_cleanup_claude.log   # net-to-1e-9 (per-shift fp64 cleanup ON)
export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES=0

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/"
LDFLAGS="-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm"
SRC=test_multishift_mixed_claude.cu

L1DIR="data_Nf6_gsq8.000000at0.200000nu01.000000nt128L1"
L2DIR="Nf6_gsq8.000000at0.200000nu01.000000nt128L2"
L4DIR="Nf6_gsq8.000000at0.200000nu01.000000nt128L4"

# 3 thermalized configs per L (spread across the stream)
L1CFGS=( "$L1DIR/F.150" "$L1DIR/F.250" "$L1DIR/F.350" )
L2CFGS=( "$L2DIR/ckpoint_lat.200" "$L2DIR/ckpoint_lat.450" "$L2DIR/ckpoint_lat.700" )
L4CFGS=( "$L4DIR/ckpoint_lat.40"  "$L4DIR/ckpoint_lat.60"  "$L4DIR/ckpoint_lat.80" )

run_L () {
  local LV=$1
  shift
  local cfgs=( "$@" )
  local bin="test_multishift_mixed_realcfg_L${LV}.x"
  echo ""
  echo "==================== BUILD LREF=${LV} Nt=128 ===================="
  if $NVCC $NVCCFLAGS $INCLUDES -DLREF=${LV} -DNTIME=128 $SRC $LDFLAGS -o "$bin" ; then
    echo ">>> BUILD OK: LREF=${LV}"
    for c in "${cfgs[@]}"; do
      echo "-------------------- RUN LREF=${LV}  config=$c --------------------"
      if [ -f "$c" ]; then ./"$bin" "$c"; else echo "!!! MISSING config $c"; fi
    done
  else
    echo "!!! BUILD FAILED: LREF=${LV}"
  fi
}

{
  echo "######## mixed-prec MULTISHIFT on REAL configs (Nf6 gsq8 massless, Nt=128)  $(date) ########"
  run_L 1 "${L1CFGS[@]}"
  run_L 2 "${L2CFGS[@]}"
  run_L 4 "${L4CFGS[@]}"
  echo ""
  echo "######## done ########"
} 2>&1 | tee "$LOG"
