#!/bin/sh
#SBATCH --account=qed3.lq2_gpu
#SBATCH --qos=test
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --time=0:30:00
#SBATCH --job-name=bench_nthreads
#SBATCH --output=/home/nmatsum/project_qed3/qed3/src/both_3d/benchmark/bench_nthreads_%j.out

BENCHDIR=/home/nmatsum/project_qed3/qed3/src/both_3d/benchmark

hostname
source /lustre2/qed3/env.sh
date
nvidia-smi

cd ${BENCHDIR}

for NT in 256 512 1024; do
  echo ""
  echo "========================================"
  echo "  NThreadsPerBlock = ${NT}"
  echo "========================================"
  time ${BENCHDIR}/bench_${NT}_claude.o
done
