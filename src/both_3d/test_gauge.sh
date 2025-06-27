#!/bin/bash -l

# https://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/
# to submit: qsub script.sh

#------ qsub options --------#
#$ -P qfe
#$ -M nmatsum@bu.edu
##### run time limit. format: hh:mm:ss; default 12 hrs
#$ -l h_rt=0:10:00


source ../../env.sh


echo "key=${key}"
dir=${key}
./a.out ${dir} plaq_ss_t_ 10000 100000000 50 ${Nt} 200 ${prefix_max} | tee "corr_${dir}.dat"

