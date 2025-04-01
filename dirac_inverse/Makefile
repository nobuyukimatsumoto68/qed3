CXX=g++
CXXFLAGS=-Wall -Wextra -std=c++23 -fuse-ld=lld -march=native -O3 -g
INCLUDES=-I"/usr/include/c++/13" -I"/usr/include/x86_64-linux-gnu/c++/13/" -I"/opt/eigen-3.4.0/" -I"/opt/gsl-2.8/"

LDFLAGS=-L"/lib/x86_64-linux-gnu/" -L"/lib/x86_64-linux-gnu/" -L"/lib64/" -L"/usr/local/lib" -lgsl -lgslcblas -lm


all:
	$(CXX) integrate.cc $(INCLUDES) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)
