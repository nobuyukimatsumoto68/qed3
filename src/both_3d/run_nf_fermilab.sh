#!/bin/sh

#SBATCH --account=qed3
#SBATCH --qos=normal
# #SBATCH --qos=test
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --time=6:00:00
# #SBATCH --time=0:10:00
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm_%x_%j.out

### # https://computing.fnal.gov/lqcd/slurm/

# print hostname of worker node
hostname
source /home/nmatsum/env.sh

# system info
date
lscpu
nvidia-smi

app=/project/qed3/hmc_fermilab_wmass_claude.o
time $app ${Nf} 1.0 ${mR} ${mI}

exit
