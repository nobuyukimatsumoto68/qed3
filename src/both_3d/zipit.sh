#!/bin/bash
#$ -P qfe

set -e

directories=$(ls -d gsq2.000000*)


# echo $directories


for directory in $directories; do
    echo $directory
    tar -czf arxv_${directory}.tar.gz $directory && rm -rf ${directory}
done
