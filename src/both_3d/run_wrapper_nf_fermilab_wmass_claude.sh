#!/bin/bash

# SRCDIR=$(pwd)
# RUNDIR=/lustre2/qed3
# Copy to rundir
# cp ${SRCDIR}/run_nf_fermilab.sh ${RUNDIR}/
# cp ${SRCDIR}/hmc_fermilab_wmass_claude.o ${RUNDIR}/
# cd ${RUNDIR}
SLURM_USER=nmatsum

affine_pairs=("0.01 0.0" "0.05 0.0" "0.0 0.01" "0.0 0.05")
qed3_pairs=("0.1 0.0" "0.2 0.0" "0.0 0.1" "0.0 0.2")

submit_job() {
    local script=$1 Nf=$2 mR=$3 mI=$4 account=$5
    local jobname="hmc_Nf${Nf}_mR${mR}_mI${mI}"
    local existing
    existing=$(squeue -u ${SLURM_USER} --name=${jobname} -h -o "%i" | tail -1)
    local dep=""
    if [ -n "${existing}" ]; then
        dep="--dependency=afterok:${existing}"
        echo "submitting (${account}) ${jobname} after job ${existing}"
    else
        echo "submitting (${account}) ${jobname}"
    fi
    sbatch --job-name=${jobname} ${dep} --export=Nf=${Nf},mR=${mR},mI=${mI} ${script}
}

for Nf in 2 4 6; do
    for pair in "${affine_pairs[@]}"; do
        mR=$(echo ${pair} | awk '{print $1}')
        mI=$(echo ${pair} | awk '{print $2}')
        submit_job run_nf_fermilab.sh ${Nf} ${mR} ${mI} affine
    done
    # for pair in "${qed3_pairs[@]}"; do
    #     mR=$(echo ${pair} | awk '{print $1}')
    #     mI=$(echo ${pair} | awk '{print $2}')
    #     submit_job run_nf_fermilab.sh ${Nf} ${mR} ${mI} qed3
    # done
done
