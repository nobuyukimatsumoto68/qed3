#!/bin/bash -l

source /projectnb/qfe/nmatsum/qed3/env.sh

# app='hmc2'
# echo $app
# export app=${app}
# qsub -N "Nf2" -v app=${app} run_nf.sh

app='hmc4'
echo $app
export app=${app}
qsub -N "Nf4" -v app=${app} run_nf.sh
