#!/bin/bash -l

# module load gcc/13.2.0
source ../../env.sh

# app='wolff.o'
app="gauge${1}"
# make ${app}

# ns=(4 6 8 10) # 1e5, omp 4
# ns=(6 8 10) # 1e5, omp 4
# len=${#ns[@]}
# for (( i=0; i<$len; i++ ));
# do
#     n=${ns[$i]}
#     dir=mult${n}
#     mkdir -p ${dir}

#     cp run.sh ${dir}
#     cp ${app} ${dir}
#     cd ${dir}
# echo $app
# export app=${app}
# if [ "$#" -eq 2 ]; then
#     qsub -N ${app} -v app=${app} -t 1-10 -j y -hold_jid ${2} run.sh
# elif [ "$#" -eq 1 ]; then
#     qsub -N ${app} -v app=${app} -t 1-10 -j y run.sh
# fi
if [ "$#" -eq 3 ]; then
    qsub -N ${app} -v app=${app} -t 1-${2} -j y -hold_jid ${3} run.sh
elif [ "$#" -eq 2 ]; then
    qsub -N ${app} -v app=${app} -t 1-${2} -j y run.sh
elif [ "$#" -eq 1 ]; then
    qsub -N ${app} -v app=${app} -j y run.sh
fi



# echo $app
# export app=${app}
# if [ "$#" -eq 2 ]; then
#     qsub -N ${app} -v app=${app} -t 1-4 -j y -hold_jid ${2} run.sh
# elif [ "$#" -eq 1 ]; then
#      qsub -N ${app} -v app=${app} -t 1-4 -j y run.sh
# fi



#     echo $n
#     cd ..
# done


# ns=(12 14) # 1e4, omp 8
# len=${#ns[@]}
# for (( i=0; i<$len; i++ ));
# do
#     n=${ns[$i]}
#     dir=mult${n}
#     mkdir -p ${dir}

#     cp run.sh ${dir}
#     cp ${app} ${dir}
#     cd ${dir}
#     qsub -N ${app}-${n} -v mult=$n -v binsize=10000 -v app=$app -j y run.sh
#     echo $n
#     cd ..
# done
