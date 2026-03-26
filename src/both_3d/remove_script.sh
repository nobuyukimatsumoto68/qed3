#!/bin/bash
set -e

# dir="./Nf4_gsq2.000000at0.200000nu01.000000nt128L1/"
dir="./data_Nf4_gsq2.000000at0.200000nu01.000000nt128L1/"

for ((i=5900; i<=9000; i++ ))
do
    # rm "${dir}ckpoint_lat.${i}"
    # rm "${dir}ckpoint_rng.${i}"
    # rm "${dir}F.${i}"
    # rm "${dir}F_corr.${i}"
done
