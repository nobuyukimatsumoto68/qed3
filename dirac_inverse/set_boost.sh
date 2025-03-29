#!/bin/bash
set -e

dir=$(pwd)/

dir_boost=${dir}/boost
if [ ! -L ${dir_boost} ] && [ ! -d ${dir_boost} ]; then
    echo "--- boost not found"
    cd ${dir}/
    # curl -O https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.bz2
    wget https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.bz2
    # git clone git@github.com:boostorg/boost.git boost_src
    tar -xf boost_1_87_0.tar.bz2
    cd boost_1_87_0
    ./bootstrap.sh
    ./b2 headers
    cd ..
    ln -s boost_1_87_0/boost boost
    echo "--- boost installed"
fi
