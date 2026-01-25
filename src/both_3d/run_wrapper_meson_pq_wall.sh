#!/bin/bash

# source /projectnb/qfe/nmatsum/qed3/env.sh

# gsq=$1
# Nf=$2
# igam=$3

app='meson_pq_wall.o'

echo $app
# if [ "$#" -eq 2 ]; then
#     echo $1 $2

make ${app}

# gsq=4.0
# Nf=2
# nu0=1.5
# nu1=${nu0}
# nhits=1
# dt=24
# ell=0
# em=0

gsq=2.0
Nf=2
nu0=3.0
nu1=3.0 # ${nu0}
nhits=1
dt=24
ell=0
em=0


for ell in 0 1 2
do
    for (( em=-ell; em<=ell; em++ ))
    do
        CUDA_VISIBLE_DEVICES=02 ./${app} ${gsq} ${Nf} ${nu0} ${nu1} ${nhits} ${dt} ${ell} ${em}
    done
done



# for gsq in 4.0
# do
#     for Nf in 2 4 6
#     do
#         # for nu0 in 0.8 1.0 1.2
#         for nu0 in 0.8 1.2
#         do
#             nu1=${nu0}
#             for ell in 0 1
#             do
#                 for (( em=-ell; em<=ell; em++ ))
#                 do
#                     echo $Nf $gsq $nu0 ${nu1} ${ell} ${em}
#                     qsub -N "mesonNf${Nf}${nu0}${ell}${em}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} -v nu0=${nu0} -v nu1=${nu1} -v ell=${ell} -v em=${em} run_meson.sh
#                 done
#             done
#         done
#     done
# done


# for gsq in 0.1 0.5 1.0 2.0
# do
#     for Nf in 2 4 6
#     do
#         echo $Nf $gsq
#         qsub -N "mesonNf${Nf}gsq${gsq}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} run_meson.sh
#     done
# done
# elif [ "$#" -eq 3 ]; then
#     echo $1 $2 $3
#     qsub -N "mesonNf${Nf}gsq${gsq}" -v app=${app} -v gsq=${gsq} -v Nf=${Nf} -hold_jid ${3} run_meson.sh
# fi



