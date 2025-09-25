CXX=g++
CXXFLAGS=-Wall -Wextra -std=c++23 -fuse-ld=lld -march=native -O3 -g
INCLUDES=-I"/usr/include/c++/13" -I"/usr/include/x86_64-linux-gnu/c++/13/" -I"/opt/eigen-3.4.0/" -I"/opt/gsl-2.8/"

qfe_dir="../include/newQFE_copy/"

QFEINCLUDES=-I"${qfe_dir}/include"

GRP_DIR="${qfe_dir}/grp/"
QFECXXFLAGS=-DGRP_DIR="\"$(GRP_DIR)\""

LDFLAGS=-L"/lib/x86_64-linux-gnu/" -L"/lib/x86_64-linux-gnu/" -L"/lib64/" -L"/usr/local/lib" -lgsl -lgslcblas -lm

all:
	$(CXX) get_geometry.cc -o geom.out $(CXXFLAGS) $(QFECXXFLAGS) $(QFEINCLUDES) $(INCLUDES) $(LDFLAGS)
	$(CXX) dual.cc -o dual.out $(INCLUDES) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)
	$(CXX) simp.cc -o simp.out $(CXXFLAGS) $(QFECXXFLAGS) $(QFEINCLUDES) $(INCLUDES) $(LDFLAGS)
