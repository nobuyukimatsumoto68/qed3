#!/bin/bash
# set -e

# Remove stuck checkpoints after the last accepted trajectory.
#
# Nf4: last accepted k=1271 (ckpoint_lat.1271 + ckpoint_rng.1269 are the restart pair).
#      Remove k=1272..1447  -> set i=1272, upper bound covers current max 1447 with margin.
#
# Nf6: last accepted k=1960 (ckpoint_lat.1960 + ckpoint_rng.1959 are the restart pair).
#      Remove k=1961..2135  -> set i=1961, upper bound covers current max 2135 with margin.
#
# To run:
#   1. Stop the running HMC process first.
#   2. Uncomment exactly ONE dir= line and set the i range accordingly (see below).
#   3. Verify with the ls lines before switching to rm.

# dir="./Nf4_gsq2.000000at0.200000nu01.000000nt128L1/"
# dir="./data_Nf4_gsq2.000000at0.200000nu01.000000nt128L1/"

# Nf4 gsq=8: remove stuck k=1272..1500  (last good: ckpoint_lat.1271 + ckpoint_rng.1269)
dir="./Nf4_gsq8.000000at0.200000nu01.000000nt128L1/"
# Nf6 gsq=8: remove stuck k=1961..2200  (last good: ckpoint_lat.1960 + ckpoint_rng.1959)
# dir="./Nf6_gsq8.000000at0.200000nu01.000000nt128L1/"
# dir="./data_Nf6_gsq8.000000at0.200000nu01.000000nt128L1/"

# Nf4: i=1272..1500   Nf6: i=1961..2200
# for ((i=1272; i<=1500; i++ ))
for ((i=1600; i<=2200; i++ ))
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
