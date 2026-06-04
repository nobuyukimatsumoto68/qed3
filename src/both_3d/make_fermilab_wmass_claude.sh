#!/bin/bash
# Build hmc_fermilab_wmass_claude.cu -> hmc_fermilab_wmass_claude.o using Makefile.fnal
# Run from anywhere; cd's to the source directory first.
source /lustre2/qed3/env.sh

# cd ${SRCDIR}
make -f Makefile.fnal hmc_fermilab_wmass_claude.o

cp hmc_fermilab_wmass_claude.o /project/qed3/hmc_fermilab_wmass_claude.o
cp hmc_fermilab_wmass_claude.o /project/affine/qed3/hmc_fermilab_wmass_claude.o
