#!/bin/bash -l
# tmp_build_distill_extend_L4_scc_claude.sh -- build-ONLY smoke test of the INCREMENTAL basis-extension binary
#   (distill_peram_extend_claude.cu: L4 Nv 48 -> 64, block width NV-NV_OLD=16) on BU SCC (sm_70 V100 + sm_80 A100/L40S).
#   Copy of tmp_build_distill_sym_L4_scc_claude.sh (untouched).  Plan: distill_peram_extend_impl_plan_claude.md.
#   Run on an SCC login node; NO submit, NO run, NO rm.
#   bash tmp_build_distill_extend_L4_scc_claude.sh 2>&1 | tee tmp_build_distill_extend_L4_scc_claude.log
set -u
cd /projectnb/qfe/nmatsum/qed3/src/production || exit 1
source /projectnb/qfe/nmatsum/qed3/env.sh
module load hdf5/1.10.10
module load gsl
NV=${NV:-64}
NV_OLD=${NV_OLD:-48}
L=4
SRC=distill_peram_extend_claude.cu
NVCCBASE="-w -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
DEFS="-DN_REFINE_CLI=${L} -DNSTACK_CLI=$(( NV - NV_OLD )) -DBASIS_SYM=1"
INCLUDES="-I./includes/ ${QED3_INC} -I/projectnb/qfe/nmatsum/opt/highfive/include -I${SCC_HDF5_INCLUDE} -I${SCC_GSL_INCLUDE}"
LDFLAGS="-L${SCC_HDF5_LIB} -L${SCC_GSL_LIB} -lhdf5 -lgsl -lgslcblas -lm"
echo "######## distill EXTEND L${L} Nv${NV_OLD}->${NV} build smoke test  $(date) ########"
ok=1
for arch in sm_70 sm_80
do
  out=distill_peram_extend_L${L}_Nv${NV_OLD}to${NV}_${arch}.out
  echo; echo "===== build $out (arch=$arch, $DEFS) ====="
  nvcc -arch="$arch" $NVCCBASE $DEFS $INCLUDES $LDFLAGS "$SRC" -o "$out" \
    && echo "  OK: $out" || { echo "  FAILED: $out"; ok=0; }
done
echo
if [ "$ok" -eq 1 ]; then echo "######## BUILD OK (both arches) ########"; else echo "######## BUILD FAILED ########"; exit 1; fi
