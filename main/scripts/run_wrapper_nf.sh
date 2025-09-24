#!/bin/bash -l

source /projectnb/qfe/nmatsum/qed3/env.sh

# app='hmc2'
# echo $app
# export app=${app}
# if [ "$#" -eq 1 ]; then
#     qsub -N "Nf2" -v app=${app} -hold_jid ${1} run_nf.sh
# else
#     qsub -N "Nf2" -v app=${app} run_nf.sh
# fi


# app='hmc4'
# echo $app
# export app=${app}
# if [ "$#" -eq 1 ]; then
#     qsub -N "Nf4" -v app=${app} -hold_jid ${1} run_nf.sh
# else
#     qsub -N "Nf4" -v app=${app} run_nf.sh
# fi



app='hmc6'
echo $app
export app=${app}
if [ "$#" -eq 1 ]; then
    qsub -N "Nf6" -v app=${app} -hold_jid ${1} run_nf.sh
else
    qsub -N "Nf6" -v app=${app} run_nf.sh
fi

