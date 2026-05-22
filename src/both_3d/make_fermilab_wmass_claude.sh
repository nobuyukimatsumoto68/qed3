#!/bin/bash
# Build hmc_fermilab_wmass_claude.cu -> hmc_fermilab_wmass_claude.o using Makefile.fnal
# Run from anywhere; cd's to the source directory first.
source /home/nmatsum/env.sh

cd ${SRCDIR}
make -f Makefile.fnal hmc_fermilab_wmass_claude.o
