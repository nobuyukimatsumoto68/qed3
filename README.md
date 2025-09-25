# qed3

All data is available in the data/ subdirectories.

## To Setup

$ source setup.sh


## Geometry data (Appendix B) for the refinement level L

This script clones a copy of an edited [copy](https://github.com/nobuyukimatsumoto68/QFE_copy/tree/main) of [newQFE](https://github.com/brower/newQFE) to generate the basic geometrical information of the refined icosahedra.

$ cd geometry

$ make

$ ./geom.out L

$ ./simp.out L

$ ./dual.out L


## Calculate eigenvalues and propagators

$ cd main

$ make

$ ./eig.o

$ ./ferm_prop.o

$ ./gauge_prop.o


## S2 check (Appendix C.3)

$ cd check_s2

$ make

$ ./ferm.out nmax

$ ./gauge.out nmax
