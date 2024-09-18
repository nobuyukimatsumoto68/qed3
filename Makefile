
CXX=g++
CXXFLAGS=-std=c++17 -g -O3

INCLUDES=-I"./newQFE/include"


GRP_DIR=$(shell pwd)/newQFE/grp/
CXXFLAGS+=-DGRP_DIR="\"$(GRP_DIR)\""


all:
	# $(CXX) get_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_geometry.cc $(INCLUDES) $(CXXFLAGS)
	# $(CXX) test_u1.cc $(INCLUDES) $(CXXFLAGS)
	$(CXX) test_dirac.cc $(INCLUDES) $(CXXFLAGS)

geometry:
	$(CXX) get_geometry.cc $(INCLUDES) $(CXXFLAGS)
