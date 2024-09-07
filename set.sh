#!/bin/bash
set -e

dir_qfe=$(pwd)/newQFE

#########################

if [ ! -L ${dir_qfe} ] && [ ! -d ${dir_qfe} ]; then
    echo "--- newQFE not found"
    git clone https://github.com/brower/newQFE
    echo "--- newQFE installed"
fi

dir_eigen=${dir_qfe}/include/Eigen
if [ ! -L ${dir_eigen} ] && [ ! -d ${dir_eigen} ]; then
    echo "--- Eigen not found"
    cd ${dir_qfe}/include/
    if [ ! -L eigen ] && [ ! -d eigen ]; then
	git clone https://gitlab.com/libeigen/eigen.git
    fi
    ln -s eigen/Eigen Eigen
    echo "--- Eigen installed"
fi

dir_boost=${dir_qfe}/include/boost
if [ ! -L ${dir_boost} ] && [ ! -d ${dir_boost} ]; then
    echo "--- boost not found"
    cd ${dir_qfe}/include/
    curl -O https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.bz2
    # wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.bz2
    # git clone git@github.com:boostorg/boost.git boost_src
    tar -xf boost_1_86_0.tar.bz2
    cd boost_1_86_0
    ./bootstrap.sh
    ./b2 headers
    cd ..
    ln -s boost_1_86_0/boost boost
    echo "--- boost installed"
fi

#########################

cd ${dir_qfe}

DBIN=${dir_qfe}/bin/

DDATA=results/q5k3/
DORBIT=results/orbit/
mkdir -p ${DDATA}
mkdir -p ${DORBIT}

FORBIT=${DORBIT}/q${q}k${k}_orbit.dat

make -C grp gen_grp_o3


cd grp
cp Makefile Makefile.mac
sed -i.bak 's/python3/${PYTHON}/g' "Makefile.mac"
export PYTHON=/opt/anaconda3/bin/python3
make gen_grp_o3 -f Makefile.mac
