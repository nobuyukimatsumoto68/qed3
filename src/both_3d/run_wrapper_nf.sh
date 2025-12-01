#!/bin/bash -l

source /projectnb/qfe/nmatsum/qed3/env.sh

gsq=0.1
app2='hmc2_01'
app4='hmc4_01'
app6='hmc6_01'

# gsq=0.5
# app2='hmc2_05'
# app4='hmc4_05'
# app6='hmc6_05'

# gsq=1.0
# app2='hmc2_10'
# app4='hmc4_10'
# app6='hmc6_10'

# gsq=2.0
# app2='hmc2_20'
# app4='hmc4_20'
# app6='hmc6_20'



export app=${app2}
echo $app
if [ "$#" -eq 3 ]; then
    qsub -N "Nf2gsq${gsq}" -v app=${app2} -hold_jid ${1} run_nf.sh
else
    qsub -N "Nf2gsq${gsq}" -v app=${app2} run_nf.sh
fi





export app=${app4}
echo $app
if [ "$#" -eq 3 ]; then
    qsub -N "Nf4gsq${gsq}" -v app=${app4} -hold_jid ${2} run_nf.sh
else
    qsub -N "Nf4gsq${gsq}" -v app=${app4} run_nf.sh
fi




export app=${app6}
echo $app
if [ "$#" -eq 3 ]; then
    qsub -N "Nf6gsq${gsq}" -v app=${app6} -hold_jid ${3} run_nf.sh
else
    qsub -N "Nf6gsq${gsq}" -v app=${app6} run_nf.sh
fi

