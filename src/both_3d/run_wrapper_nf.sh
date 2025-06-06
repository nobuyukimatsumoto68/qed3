#!/bin/bash -l

source /projectnb/qfe/nmatsum/qed3/env.sh

app='hmc.o'

echo $app
export app=${app}
qsub -N "Nf4" -v app=${app} run_nf.sh
