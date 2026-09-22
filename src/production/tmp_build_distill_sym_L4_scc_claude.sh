#!/bin/bash -l
# tmp_build_distill_sym_L4_scc_claude.sh -- build-ONLY smoke test of the SYMMETRIZED-basis distillation binary for
#   L4 on BU SCC (sm_70 V100 + sm_80 A100).  Run on an SCC login node; NO submit, NO run, NO rm.
#   bash tmp_build_distill_sym_L4_scc_claude.sh 2>&1 | tee tmp_build_distill_sym_L4_scc_claude.log
set -u
cd /projectnb/qfe/nmatsum/qed3/src/production || exit 1
source /projectnb/qfe/nmatsum/qed3/env.sh
module load hdf5/1.10.10
module load gsl
NV=${NV:-24}
L=4
SRC=distill_peram_mrhs_claude.cu
NVCCBASE="-w -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
DEFS="-DN_REFINE_CLI=${L} -DNSTACK_CLI=${NV} -DBASIS_SYM=1"
INCLUDES="-I./includes/ ${QED3_INC} -I/projectnb/qfe/nmatsum/opt/highfive/include -I${SCC_HDF5_INCLUDE} -I${SCC_GSL_INCLUDE}"
LDFLAGS="-L${SCC_HDF5_LIB} -L${SCC_GSL_LIB} -lhdf5 -lgsl -lgslcblas -lm"
echo "######## distill SYM L${L} Nv${NV} build smoke test  $(date) ########"
ok=1
for arch in sm_70 sm_80
do
  out=distill_peram_mrhs_L${L}_Nv${NV}_sym_${arch}.out
  echo; echo "===== build $out (arch=$arch, $DEFS) ====="
  nvcc -arch="$arch" $NVCCBASE $DEFS $INCLUDES $LDFLAGS "$SRC" -o "$out" \
    && echo "  OK: $out" || { echo "  FAILED: $out"; ok=0; }
done
echo
if [ "$ok" -eq 1 ]; then echo "######## BUILD OK (both arches) ########"; else echo "######## BUILD FAILED ########"; exit 1; fi
