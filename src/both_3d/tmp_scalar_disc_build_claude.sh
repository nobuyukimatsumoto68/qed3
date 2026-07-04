#!/usr/bin/env bash
# Build + free-field smoke test for the disc driver's scalar a=0 one-point loop (Chunk 3).
# Free field, m_F=0.03, disc-tblock 64 (few solves => fast smoke).  No rm anywhere.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=scalar_disc_build_claude.log
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
DISC=jj_local_ylm_scalar_disc_stoch.o

echo "### [1/4] compile scalar disc driver -> $DISC  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_scalar_disc_stoch_claude.cu -o "$DISC" || { echo "### DISC BUILD FAILED ###"; exit 1; }

echo "### [2/4] fresh disc run (m_F=0.07, disc-tblock 64) ###"
./$DISC --nhits 1 --disc-tblock 64 --mass-re 0.07 || { echo "### DISC RUN FAILED ###"; exit 1; }

echo "### [3/4] re-run --IsScalarOnly => expect SKIP (s0 already present) ###"
./$DISC --nhits 1 --disc-tblock 64 --mass-re 0.07 --IsScalarOnly || { echo "### SCALAR-ONLY RUN FAILED ###"; exit 1; }

echo "### [4/4] checks ###"
D=data_free_vmRe0.070000vmIm0.000000/corr_ylm_disc_tb64_nhits1/corr.0.h0.h5
python3 - "$D" <<'PY'
import sys
import numpy as np
import h5py
f = h5py.File(sys.argv[1], 'r')
ok = True
for k in ['h0/disc/ylm/s0/l0/m0/J/real', 'h0/disc/ylm/s0/l0/m0/J/imag',
          'h0/disc/ylm/s0/l3/m3/J/real',
          'h0/disc/ylm/s0_1mD/l0/m0/J/real', 'h0/disc/ylm/s0_1mD/l3/m3/J/imag',
          'h0/disc/ylm/s3/l1/m0/J/real', 'h0/condensate/etadag_xi/real', 'complete']:
    p = k in f; ok &= p
    print(('OK      ' if p else 'MISSING ') + k)
# s0 and s0_1mD are genuinely complex and DISTINCT (GW dressing changes the loop).
j1  = np.array(f['h0/disc/ylm/s0/l1/m0/J/real'])     + 1j*np.array(f['h0/disc/ylm/s0/l1/m0/J/imag'])
j1d = np.array(f['h0/disc/ylm/s0_1mD/l1/m0/J/real']) + 1j*np.array(f['h0/disc/ylm/s0_1mD/l1/m0/J/imag'])
print(('OK      ' if np.any(np.abs(j1.imag) > 0) else 'NOTE    ') + 's0 J has nonzero imag')
print(('OK      ' if not np.allclose(j1, j1d) else 'NOTE    ') + 's0_1mD != s0 (GW loop distinct)')
print('DISC SCALAR SMOKE:', 'PASS' if ok else 'FAIL')
PY
echo "### done  [$(date +%F_%H:%M:%S)] ###"
