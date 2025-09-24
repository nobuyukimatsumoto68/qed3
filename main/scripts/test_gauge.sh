#!/bin/bash -l

# https://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/
# to submit: qsub script.sh

#------ qsub options --------#
#$ -P qfe
#$ -M nmatsum@bu.edu
##### run time limit. format: hh:mm:ss; default 12 hrs
#$ -l h_rt=1:00:00
#$ -pe omp 8
# # # #$ -l mem_per_core=8G


source ../../env.sh


echo "key=${key}"
dir=${key}
# ./a.out ${dir} plaq_ss_t_ 10000 400000000 50 ${Nt} 2000 ${prefix_max} ${at} "corr_${dir}" | tee "corr_${dir}.dat"
# ./a.out ${dir} plaq_ss_t_ 10000 1000000000 50 ${Nt} 2000 ${prefix_max} ${at} "corr_${dir}" | tee "corr_${dir}.dat"
# ./a.out ${dir} plaq_ss_t_ 10000 100000000 50 ${Nt} 4000 ${prefix_max} ${at} "corr_${dir}" | tee "corr_${dir}.dat"

./a.out ${dir} plaq_ss_t_ 10000 10000000 50 ${Nt} 2000 ${prefix_max} ${at} "corr_${dir}" | tee "corr_${dir}.dat"

