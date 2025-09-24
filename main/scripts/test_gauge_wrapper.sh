# key="gsq0.500000at0.500000nt24L1_"
# prefix_max=1
# Nt=24
# at=0.5

# key="gsq0.500000at0.375000nt32L1_"
# prefix_max=1
# Nt=32
# at=0.375


# key="gsq0.500000at0.500000nt24L2_"
# prefix_max=8
# Nt=24
# at=0.5

# key="gsq0.500000at0.375000nt32L2_"
# prefix_max=8
# Nt=32
# at=0.375



# key="gsq0.500000at0.500000nt24L4_"
# prefix_max=16
# Nt=24
# at=0.5

key="gsq0.500000at0.375000nt32L4_"
prefix_max=16
Nt=32
at=0.375














# L=1

source ../../env.sh
g++ analyze_corr.cpp -std=c++23 -O3 # -fuse-ld=lld


# for key in ${keys[@]}
# do
export key=$key
export prefix_max=$prefix_max
export Nt=$Nt
export at=$at
echo $key
# qsub test_gauge.sh -v key=${key}
# bash test_gauge.sh -v key=${key} -v prefix_max=${prefix_max}

qsub -j y -v key=${key} -v prefix_max=${prefix_max} -v Nt=${Nt} -v at=${at} -V test_gauge.sh
# bash -v key=${key} -v prefix_max=${prefix_max} -V test_gauge.sh

# ./a.out ${dir} plaq_ss_t_ 0 10000 50 96 1000 | tee "corr_${dir}.dat"
# done




















# keys=("gsq0.050000at0.050000nt96L1_1" "gsq0.050000at0.050000nt96L1_2" "gsq0.050000at0.050000nt96L1_3" "gsq0.050000at0.050000nt96L1_4")
# keys=("gsq0.050000at0.050000nt96L2_1" "gsq0.050000at0.050000nt96L2_2" "gsq0.050000at0.050000nt96L2_3" "gsq0.050000at0.050000nt96L2_4")
# keys=("gsq0.050000at0.050000nt96L4_1" "gsq0.050000at0.050000nt96L4_2" "gsq0.050000at0.050000nt96L4_3" "gsq0.050000at0.050000nt96L4_4")



# key="gsq0.050000at0.050000nt96L1_"
# prefix_max=4



# key="gsq0.050000at0.100000nt48L2_"
# prefix_max=20
# Nt=48

# key="gsq0.050000at0.075000nt64L2_"
# prefix_max=20
# Nt=64

# key="gsq0.050000at0.050000nt96L2_"
# prefix_max=10
# Nt=96

# key="gsq0.050000at0.050000nt96L4_"
# prefix_max=20
# Nt=96

############################################



# key="gsq0.500000at0.500000nt24L1_"
# prefix_max=1
# Nt=24
# at=0.5

# key="gsq0.500000at0.375000nt32L1_"
# prefix_max=1
# Nt=32
# at=0.375

# key="gsq0.500000at0.250000nt48L1_"
# prefix_max=1
# Nt=48
# at=0.25



# key="gsq0.500000at0.500000nt24L2_"
# prefix_max=8
# Nt=24
# at=0.5

# key="gsq0.500000at0.375000nt32L2_"
# prefix_max=8
# Nt=32
# at=0.375

# key="gsq0.500000at0.250000nt48L2_"
# prefix_max=8
# Nt=48
# at=0.25



# key="gsq0.500000at0.250000nt48L4_"
# prefix_max=16/
# Nt=48
# at=0.25











# key="gsq2.000000at0.500000nt24L1_"
# prefix_max=1
# Nt=24
# at=0.5

# key="gsq2.000000at0.375000nt32L1_"
# prefix_max=1
# Nt=32
# at=0.375



# key="gsq2.000000at0.500000nt24L2_"
# prefix_max=8
# Nt=24
# at=0.5

# key="gsq2.000000at0.375000nt32L2_"
# prefix_max=8
# Nt=32
# at=0.375



# key="gsq2.000000at0.500000nt24L4_"
# prefix_max=16
# Nt=24
# at=0.5

# key="gsq2.000000at0.375000nt32L4_"
# prefix_max=16
# Nt=32
# at=0.375
