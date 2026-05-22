#!/bin/bash

gsq=8.0
Nf=4
nu0=1.0
nhits=1
dt=128
ellmax=1

run_on_gpu() {
    local gpu=$1; shift
    for nu1 in "$@"; do
        echo "=== GPU${gpu} nu1 = ${nu1} ==="
        CUDA_VISIBLE_DEVICES=${gpu} ./meson_pq_wall_v2_claude.o \
            --gsq  ${gsq}  \
            --Nf   ${Nf}   \
            --nu0  ${nu0}  \
            --nu1  ${nu1}  \
            --nhits  ${nhits}  \
            --dt     ${dt}     \
            --ellmax ${ellmax}
    done
}

# run_on_gpu 1  1.5 2.0 2.5  2>&1 | tee gpu0.log
# run_on_gpu 1  1.75 2.25 2.75  2>&1 | tee gpu1.log
run_on_gpu 0  $nu0  2>&1 | tee gpu0.log
