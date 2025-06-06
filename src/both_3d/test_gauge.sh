#!/bin/bash -l

# https://www.bu.edu/tech/support/research/system-usage/running-jobs/submitting-jobs/
# to submit: qsub script.sh

#------ qsub options --------#
#$ -P qfe
#$ -M nmatsum@bu.edu
##### run time limit. format: hh:mm:ss; default 12 hrs
#$ -l h_rt=0:30:00

source ../../env.sh

echo "${key}"
dir=${key}
./a.out ${dir} plaq_ss_t_ 50 10000 50 96 10 | tee "corr_${dir}.dat"

