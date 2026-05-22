#!/bin/sh
# # #SBATCH --account=affine
#SBATCH --account=qed3
# # #SBATCH --qos=normal
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --time=00:10:00 # HH:MM:SS

### # https://computing.fnal.gov/lqcd/slurm/


# print hostname of worker node
hostname
source /home/nmatsum/env.sh


# Run nvidia-smi in the background, logging to a file
# nvidia-smi --query-gpu=timestamp,name,utilization.gpu,memory.used --format=csv -l 30 >> gpu_usage_TH128.csv &
# Get the PID of the nvidia-smi process to kill it later if needed
# NV_PID=$!
# sleep 5
time /project/affine/nmatsum/qed3/src/both_3d/hmc_fermilab.o 12.0 6 2>&1| tee runlog_12_6_qed3.dat
# time /project/affine/nmatsum/qed3/src/both_3d/hmc_fermilab_L2.o 12.0 6 2>&1| tee runlog_L2_12_6.dat

# ... your training commands ...
# kill $NV_PID

# num_gpus=$(echo $CUDA_VISIBLE_DEVICES | awk -F',' '{print NF}')
# echo "Allocated GPUs: $num_gpus"

# num_gpus=$(nvidia-smi -L | wc -l)
# echo "Number of GPUs detected: $num_gpus"

# echo "slurm_gpus_on_node: $SLURM_GPUS_ON_NODE"
# echo "slurm_job_gpus: $SLURM_JOB_GPUS"


exit


# sacct -u nmatsum --format=JobID,JobName,AllocNodes,Elapsed,CPUTime
