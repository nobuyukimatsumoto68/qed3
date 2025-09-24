CXX=g++
CXXFLAGS=-Wall -Wextra -std=c++23 -fuse-ld=lld -march=native -O3 -g
INCLUDES=-I"/usr/include/c++/13" -I"/usr/include/x86_64-linux-gnu/c++/13/" -I"/opt/eigen-3.4.0/" -I"/opt/gsl-2.8/"

all:
	$(CXX) ferm.cc $(INCLUDES) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -o ferm.out
	$(CXX) gauge.cc $(INCLUDES) $(INCLUDES) $(LDFLAGS) $(CXXFLAGS) -o gauge.out
