#!/usr/bin/env bash
# build_sigma2_twomeson_contract_claude.sh -- build the C++/Eigen two-meson contraction (CPU only).
#   release: -O3 -DNDEBUG ; debug-ish (Eigen asserts on): -O2
set -u
cd /mnt/barracuda22/qed3/qed3/src/production || exit 1
H5=/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0
INC="-I../../qfe_mod/include -I/mnt/hdd_barracuda/opt/highfive/include -I$H5/include"
LIB="-L$H5/lib -lhdf5 -lz -ldl -lm"
SRC=sigma2_twomeson_contract_claude.cpp
g++ -std=c++20 -O3 -DNDEBUG -march=native -fopenmp $INC $SRC -o sigma2_twomeson_contract_claude.o $LIB || exit 1
g++ -std=c++20 -O2 -march=native -fopenmp $INC $SRC -o sigma2_twomeson_contract_dbg_claude.o $LIB || exit 1
echo "build OK"
