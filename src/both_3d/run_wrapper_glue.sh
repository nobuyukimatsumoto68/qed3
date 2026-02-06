#!/bin/bash

source /projectnb/qfe/nmatsum/qed3/env.sh

# gsq=$1
# Nf=$2
# igam=$3

app='glue2.o'

echo $app
# if [ "$#" -eq 2 ]; then
#     echo $1 $2


for gsq in 2.8
do
    for Nf in 2 4 6 # 4 6
    do
        for nu0 in 1.0
        do
            echo $Nf $gsq $nu0
            ./${app} ${gsq} ${Nf} ${nu0}
            # qsub -N "glueNf${Nf}gsq${gsq}${nu0}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} -v nu0=${nu0} run_glue.sh
        done
    done
done
# elif [ "$#" -eq 3 ]; then
#     echo $1 $2 $3
#     qsub -N "mesonNf${Nf}gsq${gsq}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} -hold_jid ${3} run_meson.sh
# fi



