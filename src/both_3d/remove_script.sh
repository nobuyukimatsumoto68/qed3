#!/bin/bash
# set -e

# dir="./Nf4_gsq2.000000at0.200000nu01.000000nt128L1/"
# dir="./data_Nf4_gsq2.000000at0.200000nu01.000000nt128L1/"

dir="./Nf4_gsq8.000000at0.200000nu01.000000nt128L1/"
# dir="./Nf6_gsq8.000000at0.200000nu01.000000nt128L1/"
# dir="./data_Nf6_gsq8.000000at0.200000nu01.000000nt128L1/"

for ((i=1200; i<=2600; i++ ))
do
    rm -v "${dir}ckpoint_lat.${i}"
    rm -v "${dir}ckpoint_rng.${i}"
    # ls "${dir}ckpoint_lat.${i}"
    # ls "${dir}ckpoint_rng.${i}"
    # rm "${dir}F.${i}"
    # rm "${dir}F_corr.${i}"
    # ls "${dir}F.${i}"
    # ls "${dir}F_corr.${i}"
done
