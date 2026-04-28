#!/bin/sh
#SBATCH --account=affine
# # #SBATCH --qos=normal
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1

### # https://computing.fnal.gov/lqcd/slurm/


# print hostname of worker node
hostname
source /home/nmatsum/env.sh
# sleep 5
time /project/affine/nmatsum/qed3/src/both_3d/meson_pq_wall_v2.o --gsq 12.0 --Nf 6 --nu0 1.0 --nu1 1.0 --nhits 1 --dt 192 --ellmax 2 2>&1| tee log_meson_12_6.dat
# time /project/affine/nmatsum/qed3/src/both_3d/meson_pq_wall_v2_L2.o --gsq 12.0 --Nf 6 --nu0 1.0 --nu1 1.0 --nhits 1 --dt 192 --ellmax 2 2>&1| tee log_meson_L2_12_6.dat
# time /project/affine/nmatsum/qed3/src/both_3d/meson_pq_wall_v2_L2.o 12.0 6 2>&1| tee log_meson_L2_12_6.dat
exit


# sacct -u nmatsum --format=JobID,JobName,AllocNodes,Elapsed,CPUTime
