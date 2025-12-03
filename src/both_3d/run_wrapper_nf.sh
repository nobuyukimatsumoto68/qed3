#!/bin/bash -l

source /projectnb/qfe/nmatsum/qed3/env.sh

gsq=$1
Nf=$2

app='hmc.o'

echo $app
if [ "$#" -eq 2 ]; then
    echo $1 $2
    qsub -N "Nf${Nf}gsq${gsq}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} run_nf.sh
elif [ "$#" -eq 3 ]; then
    echo $1 $2 $3
    qsub -N "Nf${Nf}gsq${gsq}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} -hold_jid ${3} run_nf.sh
fi



