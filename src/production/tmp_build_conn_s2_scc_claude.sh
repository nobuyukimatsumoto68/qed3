#!/bin/bash -l
# tmp_build_conn_s2_scc_claude.sh  (_claude handoff, 2026-08-10)
# Build-ONLY smoke test for the SCC L4 conn stride-2 binary (both arches). NM runs this on an SCC login
# node; Claude reads back build_conn_s2_L4_*_claude.log. NO submit, NO run, NO rm.
#   bash tmp_build_conn_s2_scc_claude.sh 2>&1 | tee tmp_build_conn_s2_scc_claude.log
set -u

cd /projectnb/qfe/nmatsum/qed3/src/production || exit 1
source /projectnb/qfe/nmatsum/qed3/env.sh
module load hdf5/1.10.10
module load gsl

echo "######## conn-s2 L4 build smoke test  $(date) ########"
echo "# QED3_INC=$QED3_INC"
echo "# SCC_HDF5_INCLUDE=$SCC_HDF5_INCLUDE"
echo "# SCC_GSL_INCLUDE=$SCC_GSL_INCLUDE"
echo "# HighFive=/projectnb/qfe/nmatsum/opt/highfive/include"

SRC=jj_local_ylm_scalar_conn_stoch_claude.cu
NVCCBASE="-g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES="-I./includes/ ${QED3_INC} -I/projectnb/qfe/nmatsum/opt/highfive/include -I${SCC_HDF5_INCLUDE} -I${SCC_GSL_INCLUDE}"
LDFLAGS="-L${SCC_HDF5_LIB} -L${SCC_GSL_LIB} -lhdf5 -lgsl -lgslcblas -lm"

ok=1
for arch in sm_70 sm_80
do
  out=jj_conn_s2_L4_${arch}.out
  echo; echo "===== build $out (arch=$arch, -DN_REFINE_CLI=4) ====="
  nvcc -arch="$arch" $NVCCBASE -DN_REFINE_CLI=4 $INCLUDES $LDFLAGS "$SRC" -o "$out" \
    && echo "  OK: $out" || { echo "  FAILED: $out"; ok=0; }
done

echo
if [ "$ok" -eq 1 ]
then
  echo "######## BUILD OK (both arches) -- next: bash run_wrapper_conn_s2_scc_claude.sh ########"
else
  echo "######## BUILD FAILED -- see errors above ########"
fi
