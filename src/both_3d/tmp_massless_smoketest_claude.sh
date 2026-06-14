#!/usr/bin/env bash
# MASSLESS (non-parity) REGRESSION smoke test for the parity refactor (audit I1-I4).  The refactor only
# touched the if(parity) blocks, so the massless/m_F path must be BYTE-UNCHANGED.  Free U=1, mass=0,
# nhits=1, n-t0=2, spin x td=2 (4 patterns), GPU 0.  Distinct binary jj_corr_dilute_msltest.o.
#
# Two checks:
#  (A) key layout: massless writes h0/{tp,sp}/Vmm (= conj(Vpp)), h0/s{c}/Vmm, h0/disc/{tp,sp}/J, condensate;
#      and must NOT have the parity-only h0/disc/{tp,sp}/Jtil, nor any h0/t0_<b>/.../Vmm nesting.
#  (B) regression: compare the new hit-0 file to the EXISTING massless reference (same RNG seed -> should
#      match) corr_dil_nt02_nhits4_s1td2/corr.0.h0.h5 (built pre-refactor by tmp_condensate_deter).
#      max |diff| over all h0/ datasets should be ~roundoff (<1e-10) if the non-parity path is unchanged.
#
# Usage: bash tmp_massless_smoketest_claude.sh [GPU]   (default 0).  Read back: massless_smoketest_claude.log
# NOTE: if a COMPLETE corr.0.h0.h5 already exists in the nhits1 output dir it is resume-skipped (rm by hand).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1
export CUDA_VISIBLE_DEVICES="${1:-0}"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
LOG=massless_smoketest_claude.log

NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

NEW=data_free_vmRe0.000000vmIm0.000000/corr_dil_nt02_nhits1_s1td2/corr.0.h0.h5
REF=data_free_vmRe0.000000vmIm0.000000/corr_dil_nt02_nhits4_s1td2/corr.0.h0.h5   # pre-refactor reference (hit 0)

echo "### compile -> jj_corr_dilute_msltest.o  [$(date +%F_%H:%M:%S)] ###" | tee "$LOG"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o jj_corr_dilute_msltest.o 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### BUILD FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### run FREE massless nhits=1 n-t0=2 spin+td2 (GPU $CUDA_VISIBLE_DEVICES)  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
./jj_corr_dilute_msltest.o --nhits 1 --n-t0 2 --spin-dilution --time-dilution 2 2>&1 | tee -a "$LOG"
st=${PIPESTATUS[0]}; if [ "$st" -ne 0 ]; then echo "### RUN FAILED ($st) ###" | tee -a "$LOG"; exit 1; fi

echo "### (A) key-layout + (B) regression vs $REF  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
python3 - "$NEW" "$REF" <<'PY' 2>&1 | tee -a "$LOG"
import sys, h5py, numpy as np
new, ref = sys.argv[1], sys.argv[2]
need   = ['h0/tp/Vpp','h0/tp/Vmm','h0/sp/Vpp','h0/sp/Vmm','h0/s1/Vmm','h0/s2/Vmm','h0/s3/Vmm',
          'h0/axial/tp/Apm','h0/disc/tp/J','h0/disc/sp/J',
          'h0/condensate/etadag_xi','h0/condensate/xidag_eta','h0/condensate/xidag_1mDdag_eta']
absent = ['h0/disc/tp/Jtil','h0/disc/sp/Jtil','h0/t0_0/tp/Vmm','h0/t0_1/tp/Vmm']  # parity-only / old nesting
ok = True
with h5py.File(new,'r') as f:
    print('complete flag :', 'complete' in f); ok = ok and ('complete' in f)
    for k in need:   p=(k+'/real') in f; print('  present', k, ':', p); ok = ok and p
    for k in absent: a=(k+'/real') not in f; print('  absent ', k, ':', a); ok = ok and a
print('(A) key layout:', 'PASS' if ok else 'FAIL')

# (B) regression: max |diff| over all leaf datasets under h0/
import os
regr = None
if os.path.exists(ref):
    def leaves(g,p=''):
        out={}
        for k in g:
            it=g[k]; np_=p+'/'+k
            if isinstance(it,h5py.Group): out.update(leaves(it,np_))
            else: out[np_]=it[()]
        return out
    with h5py.File(new) as a, h5py.File(ref) as b:
        la, lb = leaves(a['h0']), leaves(b['h0'])
        common = sorted(set(la) & set(lb))
        mx = 0.0; worst=''
        for k in common:
            d = float(np.max(np.abs(np.asarray(la[k]).ravel() - np.asarray(lb[k]).ravel()))) if np.asarray(la[k]).size else 0.0
            if d>mx: mx=d; worst=k
        regr = mx
        print('(B) regression vs pre-refactor hit-0: max|diff|=%.3e (worst %s) over %d datasets'%(mx,worst,len(common)))
        print('    -> non-parity path UNCHANGED' if mx<1e-10 else '    -> DIFFERS (>1e-10): investigate')
else:
    print('(B) regression: reference %s not found -- skipped'%ref)
sys.exit(0 if (ok and (regr is None or regr<1e-10)) else 2)
PY
echo "### massless smoke test done (status ${PIPESTATUS[0]})  [$(date +%F_%H:%M:%S)] ###" | tee -a "$LOG"
