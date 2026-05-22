#!/bin/sh
#SBATCH --account=affine.lq2_gpu
#SBATCH --qos=normal
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --time=6:00:00
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

# GPU monitoring in background
# nvidia-smi --query-gpu=timestamp,name,utilization.gpu,memory.used --format=csv -l 30 > gpu_usage_${Nf}_${mR}_${mI}.csv &
# NV_PID=$!

time /project/affine/nmatsum/qed3/src/both_3d/hmc_fermilab_wmass_claude.o 8.0 ${Nf} 1.0 ${mR} ${mI}

# kill $NV_PID

exit


# sacct -u nmatsum --format=JobID,JobName,AllocNodes,Elapsed,CPUTime
