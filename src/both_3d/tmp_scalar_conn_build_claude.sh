#!/usr/bin/env bash
# Build + free-field smoke test for the scalar connected ylm driver (FURNISHED: h0/scalar + h0/scalar_fs).
# Phase A: default mode (vector+axial+scalar in one fresh file), m_F=0.07.
# Phase B: --IsScalarOnly APPEND path -- build the ORIGINAL conn driver to make a vector+axial-only file
#          (m_F=0.09), then append h0/scalar + h0/scalar_fs and verify the vector keys survive + scalar added.
# Fresh masses (0.07/0.09) => distinct dirs, no collision with the earlier identity-version smoke.  No rm.
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d

LOG=scalar_conn_build_claude.log
: > "$LOG"
exec > >(tee -a "$LOG") 2>&1

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'
SCAL=jj_local_ylm_scalar_conn_stoch.o
CONN=jj_local_ylm_conn_stoch.o

echo "### [1/5] compile scalar driver -> $SCAL  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_scalar_conn_stoch_claude.cu -o "$SCAL" || { echo "### SCALAR BUILD FAILED ###"; exit 1; }

echo "### [2/5] Phase A: default free run (m_F=0.07, spin dilution, 1 hit) ###"
./$SCAL --nhits 1 --t0 0 --spin-dilution --mass-re 0.07 || { echo "### RUN A FAILED ###"; exit 1; }

echo "### [3/5] compile ORIGINAL conn driver -> $CONN (for the append fixture) ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o "$CONN" || { echo "### CONN BUILD FAILED ###"; exit 1; }

echo "### [4/5] Phase B: make vector+axial-only file (m_F=0.09) then APPEND scalar+scalar_fs ###"
./$CONN --nhits 1 --t0 0 --spin-dilution --mass-re 0.09 || { echo "### CONN RUN FAILED ###"; exit 1; }
./$SCAL --nhits 1 --t0 0 --spin-dilution --mass-re 0.09 --IsScalarOnly || { echo "### APPEND RUN FAILED ###"; exit 1; }

echo "### [5/5] checks ###"
A=data_free_vmRe0.070000vmIm0.000000/corr_ylm_conn_t00_nhits1_s1/corr.0.h0.h5
B=data_free_vmRe0.090000vmIm0.000000/corr_ylm_conn_t00_nhits1_s1/corr.0.h0.h5
python3 - "$A" "$B" <<'PY'
import sys
import numpy as np
import h5py

def has(f, k):
    return k in f

ok = True

# Phase A: fresh full file has scalar (PS) + scalar_fs (FS V--) + vector + complete; Vmm == conj(Vpp).
fa = h5py.File(sys.argv[1], 'r')
for k in ['h0/scalar/l0/m0/Vpp/real', 'h0/scalar/l3/m3/Vpp/real',
          'h0/scalar_fs/l1/m0/Vmm/real', 'h0/scalar_fs/l3/m-3/Vmm/real',
          'h0/ylm/s3/l1/m0/Vpp/real', 'complete']:
    p = has(fa, k); ok &= p
    print(('A OK      ' if p else 'A MISSING ') + k)
vpp_i = np.array(fa['h0/scalar/l1/m0/Vpp/imag'])
vmm_i = np.array(fa['h0/scalar/l1/m0/Vmm/imag'])
c = np.allclose(vmm_i, -vpp_i); ok &= c
print(('A OK      ' if c else 'A FAIL    ') + 'Vmm == conj(Vpp)')
# scalar_fs (V--^FS) is genuinely GW-dressed => NOT equal to conj(Vpp) (distinguishes it from the PS Vmm).
fs_r = np.array(fa['h0/scalar_fs/l1/m0/Vmm/real'])
vpp_r = np.array(fa['h0/scalar/l1/m0/Vpp/real'])
distinct = not np.allclose(fs_r, vpp_r)
print(('A OK      ' if distinct else 'A NOTE    ') + 'scalar_fs V-- != Vpp (GW dressing present)')

# Phase B: appended file KEEPS vector+axial+complete AND gains scalar + scalar_fs.
fb = h5py.File(sys.argv[2], 'r')
for k in ['h0/ylm/s3/l1/m0/Vpp/real', 'h0/ylm_axial/s1/l2/m0/Vpp/real', 'complete',
          'h0/scalar/l0/m0/Vpp/real', 'h0/scalar_fs/l2/m0/Vmm/real']:
    p = has(fb, k); ok &= p
    print(('B OK      ' if p else 'B MISSING ') + k)

print('SCALAR SMOKE:', 'PASS' if ok else 'FAIL')
PY
echo "### done  [$(date +%F_%H:%M:%S)] ###"
