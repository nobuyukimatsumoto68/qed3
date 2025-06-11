# keys=("gsq0.050000at0.050000nt96L1_1" "gsq0.050000at0.050000nt96L1_2" "gsq0.050000at0.050000nt96L1_3" "gsq0.050000at0.050000nt96L1_4")
# keys=("gsq0.050000at0.050000nt96L2_1" "gsq0.050000at0.050000nt96L2_2" "gsq0.050000at0.050000nt96L2_3" "gsq0.050000at0.050000nt96L2_4")
# keys=("gsq0.050000at0.050000nt96L4_1" "gsq0.050000at0.050000nt96L4_2" "gsq0.050000at0.050000nt96L4_3" "gsq0.050000at0.050000nt96L4_4")


key="gsq0.050000at0.050000nt96L1_"
# key="gsq0.050000at0.050000nt96L2_"
# key="gsq0.050000at0.050000nt96L4_"
prefix_max=4

L=1

source ../../env.sh
g++ analyze_corr.cpp -std=c++23 -O3 -fuse-ld=lld


# for key in ${keys[@]}
# do
export key=$key
export prefix_max=$prefix_max
echo $key
# qsub test_gauge.sh -v key=${key}
bash test_gauge.sh -v key=${key} -v prefix_max=${prefix_max}
# ./a.out ${dir} plaq_ss_t_ 0 10000 50 96 1000 | tee "corr_${dir}.dat"
# done
