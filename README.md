# qed3


## To Setup

$ source setup.sh


## Geometry data (Appendix B) for the refinement level L

$ cd geometry

$ make

$ ./geom.out L

$ ./simp.out L

$ ./dual.out L


## Eigenvalues and propagators


## S2 check (Appendix C.3)

$ cd check_s2

$ make

$ ./ferm.out nmax

$ ./gauge.out nmax
