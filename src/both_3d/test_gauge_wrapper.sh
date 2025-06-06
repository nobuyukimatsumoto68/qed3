keys=("gsq0.050000at0.050000nt96L1_1" "gsq0.050000at0.050000nt96L1_2" "gsq0.050000at0.050000nt96L1_3" "gsq0.050000at0.050000nt96L1_4")

L=1
suffix=1

source ../../env.sh
g++ analyze_corr.cpp -std=c++23 -O3 -fuse-ld=lld


for key in ${keys[@]}
do
    export key=$key
    echo $key
    qsub test_gauge.sh -v key=${key}
    # ./a.out ${dir} plaq_ss_t_ 0 10000 50 96 1000 | tee "corr_${dir}.dat"
done
