#!/bin/bash -l

# https://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/
# to submit: qsub script.sh

#------ qsub options --------#
#$ -P qfe
#$ -M nmatsum@bu.edu
##### run time limit. format: hh:mm:ss; default 12 hrs
#$ -l h_rt=0:30:00

module load gcc/13.2.0
g++ analyze_corr.cpp -std=c++17 -O3

# key=beta50.000000at0.100000
# key=beta20.000000at0.100000
# key=beta5.000000at0.100000
# key=beta5.000000at0.200000

# key=beta20.000000at0.025000
# key=beta20.000000at0.025000
# key=beta50.000000at0.050000
# key=beta20.000000at0.050000
# key=beta50.000000at0.100000

# key=beta50.000000at0.100000
# key=beta50.000000at0.050000
# key=beta20.000000at0.050000
key=beta20.000000at0.025000

# drwxr-sr-x 2 nmatsum qfe   2097152 May 12 10:28 beta50.000000at0.050000nt96L2_short
# drwxr-sr-x 2 nmatsum qfe   2097152 May 12 10:28 beta20.000000at0.025000nt96L2_short
# drwxr-sr-x 2 nmatsum qfe   2097152 May 12 10:28 beta20.000000at0.050000nt96L2_short
# drwxr-sr-x 2 nmatsum qfe    524288 May 12 10:28 beta20.000000at0.050000nt96L4_short
# drwxr-sr-x 2 nmatsum qfe    524288 May 12 10:27 beta20.000000at0.025000nt96L4_short
# drwxr-sr-x 2 nmatsum qfe    262144 May 12 10:27 beta50.000000at0.050000nt96L4_short
L=1

dir=${key}nt96L${L}_short
./a.out ${dir} plaq_ss_t_ 0 100001000 20 96 1000 | tee "corr_${dir}.dat"

# dir=beta20.000000at0.100000nt96L2
# ./a.out ${dir} plaq_ss_t_ 0 10001000 20 96 100

# dir=beta20.000000at0.100000nt96L4
# ./a.out ${dir} plaq_ss_t_ 0 10001000 20 96 100
