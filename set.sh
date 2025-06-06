#!/bin/bash
set -e

dir_qfe=$(pwd)/../qfe_mod

#########################

if ! [ -x "$(command -v curl)" ]; then
    echo 'Error: curl is not installed. Please install with, e.g., brew install curl'
    exit 1
fi

#########################

if [ ! -L ${dir_qfe} ] && [ ! -d ${dir_qfe} ]; then
    echo "--- qfe_mod not found"
    git clone https://github.com/nobuyukimatsumoto68/qfe_mod.git
    echo "---  qfe_mod installed"
fi

dir_eigen=${dir_qfe}/include/Eigen
if [ ! -L ${dir_eigen} ] && [ ! -d ${dir_eigen} ]; then
    echo "--- Eigen not found"
    cd ${dir_qfe}/include/
    if [ ! -L eigen ] && [ ! -d eigen ]; then
	git clone https://gitlab.com/libeigen/eigen.git eigen_src
    fi
    ln -s eigen_src/Eigen Eigen
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

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    make -C grp gen_grp_o3
    # make -C grp gen_grp_o3_V2
elif [[ "$OSTYPE" == "darwin"* ]]; then
    cd ${dir_qfe}/grp
    cp Makefile Makefile.mac
    sed -i.bak 's/python3/${PYTHON}/g' "Makefile.mac"

    PATHTOCONDA=$(conda info | grep "base environment" | cut -f 2 -d : | cut -f 1 -d "(" | sed "s/ //g")
    if [[ -x "${PATHTOCONDA}/bin/python3" ]]; then
	export PYTHON=${PATHTOCONDA}/bin/python3
    elif [[ -x "$(command -v python3)" ]]; then
	export PYTHON=$(command -v python3)
    elif [[ -x "$(command -v python)" ]]; then
	export PYTHON=$(command -v python)
    else
	echo 'Error: python is not installed. Please install with, e.g., brew install python3'
	exit 1
    fi

    cd ${dir_qfe}
    make gen_grp_o3 -C grp -f Makefile.mac
else
    echo "unknown OS"
    exit 1
fi
