#!/bin/bash -l

source /projectnb/qfe/nmatsum/qed3/env.sh

gsq=$1
Nf=$2
nu0=$3

app='hmc.o'

my_array=()

echo $app
# if [ "$#" -eq 3 ]; then
#     echo $1 $2 $3
for gsq in 2.4
do
    for Nf in 6
    # for Nf in 4
    do
        for nu0 in 1.0
        # for nu0 in 1.0 1.2
        # for nu0 in 1.0
        do
            echo $Nf $gsq $nu0
            CUDA_VISIBLE_DEVICES=0 ./${app} ${gsq} ${Nf} ${nu0}
            # qsub -N "Nf${Nf}gsq${gsq}nu${nu0}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} -v nu0=${nu0} run_nf.sh
            # qsub -N "Nf${Nf}gsq${gsq}nu0${nu0}" -V run_nf.sh
            # elif [ "$#" -eq 4 ]; then
            # echo $1 $2 $3 $4
            # qsub -N "Nf${Nf}gsq${gsq}nu0${nu0}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} -v nu0=${nu0} -hold_jid ${4} run_nf.sh
        done
    done
done


