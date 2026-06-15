#!/usr/bin/env python3
# Build the three per-mass validation notebooks (m0, mF, mP), each with vector + axial + condensate.
# Single source of truth -> guarantees parallel structure.  Conventions (user, 2026-06-15):
#   - Nt/2-shift translation average everywhere (half_shift_avg); NO four-fold averaging.
#   - LOG-scale correlator plots, with the two signs in SEPARATE cells: each Re/Im part gives a
#     '+X on log' cell and a '-X on log' cell (the sign branch with data shows the decay; the other is empty).
#   - independent region 0<t<Nt/2.  compare diluted vs the determ ground truth only (no CFT-shape overlay).
#   - Vmm (V--) cells: mP notebook ONLY.  Local cells plot the Vpp channel only.
# Output: jj_validate_{m0,mF,mP}_claude.ipynb .  Re-run to regenerate.  Does NOT delete anything.
import json, uuid

# ---------------------------------------------------------------- per-mass config cell (token __MASS__)
CONFIG = r'''# ----- config (__MASS__) -----
import numpy as np, h5py, glob, os
import matplotlib.pyplot as plt

dir1 = '/mnt/barracuda22/qed3/qed3/src/both_3d/'
MASS = '__MASS__'                       # 'massless' | 'mF' | 'mP'  (this notebook is dedicated to one mass)
_ESN = {'massless':'data_free_vmRe0.000000vmIm0.000000/',
        'mF':      'data_free_vmRe0.100000vmIm0.000000/',
        'mP':      'data_free_vmRe0.000000vmIm0.100000/'}
ESN  = dir1 + _ESN[MASS]
at_cft = 0.2; Nt_cft = 128
DIL   = 'corr_dil_nt02_nhits140_s1td2'  # 140-hit validation run (spin x e/o-time dilution = 4 patterns/hit)
VOL   = 'corr_nt02_nhits64'             # volume run (only exists for m0; variance cell auto-skips otherwise)
LOC   = 'corr_deter_local_L1'           # local vector determ
LOCAX = 'corr_deter_local_axial_L1'     # local axial determ
EXS1A = 'corr_deter_exact1_axial_L1'    # exact-K axial determ (single insertion)

# vector exact-K determ: m0 has the full exactsum (x1); mF/mP have only exact1 (exactsum = 4pi * exact1).
if MASS == 'massless':
    EXS = 'corr_deter_exactsum_L1'; EXS_SCALE = 1.0
else:
    EXS = 'corr_deter_exact1_L1';   EXS_SCALE = 4.0*np.pi
CFAC = 1.0/(4.0*np.pi)                  # axial: dilute = area-weighted SUM w/ 1/4pi vs exact1 single insertion
parity = (MASS == 'mP')                 # parity (m_P): disc Vmm uses Jtil; else Vmm = conj(Vpp)

_MLAB = {'massless':'m0','mF':'m0.1','mP':'m0.1i'}; MLAB = _MLAB[MASS]   # savefig filename tag
_MTIT = {'massless':'0', 'mF':'0.1', 'mP':'0.1i'};  MTIT = _MTIT[MASS]   # title tag ($m=...$)
FIGDIR = 'figs'; os.makedirs(FIGDIR, exist_ok=True)
print('MASS=%s  ESN=%s  EXS=%s (x%.4f)' % (MASS, _ESN[MASS], EXS, EXS_SCALE))
'''

HELPERS = r'''# ----- shared helpers -----
# Folding: ONLY the master-field Nt/2-shift TRANSLATION average (period Nt/2); NO four-fold reflection fold.
# Correlator plots are LOG-scale with the two signs in separate cells (+X and -X), 0<t<Nt/2.
def dil_files(dirn):
    out=[]
    for f in sorted(glob.glob(ESN+dirn+'/corr.0.h*.h5')):
        with h5py.File(f,'r') as h:
            if 'complete' in h: out.append(f)
    return out
def dil_hits(dirn, key):                 # key e.g. 'tp/Vpp','s3/Vpp','axial/tp/Apm' -> (nhits,Nt) complex
    out=[]
    for f in dil_files(dirn):
        with h5py.File(f,'r') as h:
            out.append(h['h0/'+key+'/real'][()] + 1j*h['h0/'+key+'/imag'][()])
    assert out, 'no complete diluted hits in %s' % dirn
    return np.array(out)
def jk_cplx(arr):                        # complex (nhits,M) -> (mean, err_re, err_im)
    H=len(arr); s=arr.sum(0)
    samp=np.array([(s-arr[i])/(H-1) for i in range(H)])
    est=samp.mean(0)
    vre=(H-1)*np.mean((samp.real-est.real)**2,0); vim=(H-1)*np.mean((samp.imag-est.imag)**2,0)
    return est, np.sqrt(vre), np.sqrt(vim)
def half_shift_avg(C):                   # Nt/2-shift translation average (mean kept, stoch noise ~/sqrt2)
    return 0.5*(C + np.roll(C, C.shape[-1]//2, axis=-1))
def mirror(G):                           # single-origin G -> master-field G + roll(G, Nt/2)
    return G + np.roll(G, len(G)//2)
def det_raw(corrdir, sub):               # single-origin determ, sub e.g. 'tp/Vpp','tp/Apm','s3/Vpp' -> complex (Nt,)
    with h5py.File(ESN+corrdir+'/corr.0.h5','r') as f:
        b='h0/t0_0/'+sub
        return f[b+'/real'][()] + 1j*f[b+'/imag'][()]
def jk_ratio_re(num, den):               # jackknife ratio of Re(half-shift) hit-sums -> (mean,err)
    n=num.real; d=den.real; H=len(n); sn=n.sum(0); sd=d.sum(0)
    samp=np.array([(sn-n[i])/(sd-d[i]) for i in range(H)])
    est=samp.mean(0); var=(H-1)*np.mean((samp-est)**2,0)
    return est, np.sqrt(var)

# vector exact (Vpp/Vmm): diluted (jk) + determ in the diluted (area-weighted-summed) normalization.
def prep_v(proj, pm):
    m,ere,eim = jk_cplx(half_shift_avg(dil_hits(DIL, proj+'/'+pm)))
    d = half_shift_avg(mirror(EXS_SCALE*det_raw(EXS, proj+'/'+pm)))
    return m,ere,eim,d
# vector local s{c} (Vpp channel): diluted (jk) + determ (direct, same summed norm).
def prep_loc(c):
    m,ere,eim = jk_cplx(half_shift_avg(dil_hits(DIL, 's%d/Vpp'%c)))
    d = half_shift_avg(mirror(det_raw(LOC, 's%d/Vpp'%c)))
    return m,ere,eim,d
# vector local G_s = (s1+s2)/2 (Vpp): diluted (jk) + determ.
def prep_locsp():
    A = 0.5*(half_shift_avg(dil_hits(DIL,'s1/Vpp'))+half_shift_avg(dil_hits(DIL,'s2/Vpp')))
    m,ere,eim = jk_cplx(A)
    d = 0.5*(half_shift_avg(mirror(det_raw(LOC,'s1/Vpp')))+half_shift_avg(mirror(det_raw(LOC,'s2/Vpp'))))
    return m,ere,eim,d
# axial exact tp/sp: diluted (jk) + determ exact1 + REAL matching const c (=1/4pi: summed vs single insertion).
def prep_ax(proj):
    m,ere,eim = jk_cplx(half_shift_avg(dil_hits(DIL, 'axial/%s/Apm'%proj)))
    d = half_shift_avg(mirror(det_raw(EXS1A, '%s/Apm'%proj)))
    return m,ere,eim,d,CFAC
# diluted local axial s{c} -> (nhits,Nt) complex, or None if absent
def dil_ax_loc(c):
    out=[]
    for f in dil_files(DIL):
        with h5py.File(f,'r') as h:
            k='h0/axial/s%d/Apm'%c
            if (k+'/real') not in h: return None
            out.append(h[k+'/real'][()]+1j*h[k+'/imag'][()])
    return np.array(out) if out else None
# local axial s{c}: diluted (jk) + determ (direct), or None
def prep_axloc(c):
    A = dil_ax_loc(c)
    if A is None: return None
    m,ere,eim = jk_cplx(half_shift_avg(A))
    d = half_shift_avg(mirror(det_raw(LOCAX, 's%d/Apm'%c)))
    return m,ere,eim,d
print('diluted hits in %s: %d' % (DIL, len(dil_files(DIL))))
'''

# ---- one signed-log cell per Re/Im part: BOTH sign branches in the same plot (prep sets m,ere,eim,d[,c]).
# On a log axis the positive points are drawn as filled markers, |negative| as open markers; the determ
# curve likewise splits into a '+' (solid) and '-' (dotted) branch.  Plot code is inline per the house style.
TEMPLATE = (
"# @TITLE@ @P@ on signed-log (both sign branches; diluted vs determ).  0<t<Nt/2.\n"
"@PREP@"
"v = @VAL@; e = @ERR@; dv = @DVAL@\n"
"pos = v[dt] > 0; neg = v[dt] < 0; dpos = dv[dt] > 0; dneg = dv[dt] < 0\n"
"plt.figure()\n"
"plt.errorbar(t[pos], v[dt][pos], yerr=e[dt][pos], fmt='o', ms=3, capsize=2, label='diluted $+$')\n"
"plt.errorbar(t[neg], -v[dt][neg], yerr=e[dt][neg], fmt='v', ms=3, capsize=2, mfc='none', label='diluted $-$')\n"
"plt.plot(t[dpos], dv[dt][dpos], 's--', ms=3, label='determ $+$')\n"
"plt.plot(t[dneg], -dv[dt][dneg], 'x:', ms=4, label='determ $-$')\n"
"plt.yscale('log'); plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$|\\mathrm{@P@}\\,@YLAB@|$')\n"
"plt.title('@TITLE@ @P@ (signed log), $m='+MTIT+'$'); plt.legend(fontsize=8)\n"
"plt.savefig(FIGDIR+'/@FN@_'+MLAB+'.pdf', bbox_inches='tight')\n"
)

def cell(prep, part, ylab, title, fname, cfac=False):
    P   = 'Re' if part=='re' else 'Im'
    suf = '.real' if part=='re' else '.imag'
    val = ('(c*m)' if cfac else 'm') + suf
    err = ('c*' if cfac else '') + ('ere' if part=='re' else 'eim')
    dval= 'd' + suf
    T = TEMPLATE
    T = T.replace('@PREP@', prep).replace('@VAL@', val).replace('@ERR@', err).replace('@DVAL@', dval)
    T = T.replace('@YLAB@', ylab).replace('@TITLE@', title).replace('@FN@', fname+'_'+part).replace('@P@', P)
    return T

# ---- guarded signed-log cell for the local axial (prep_axloc may return None) ----
def axloc_cell(c, part):
    P   = 'Re' if part=='re' else 'Im'
    suf = '.real' if part=='re' else '.imag'
    err = 'ere' if part=='re' else 'eim'
    return ("# axial local s%d %s on signed-log (both sign branches; diluted vs determ).  0<t<Nt/2.\n" % (c, P)
            + "dt = np.arange(1, Nt_cft//2); t = dt*at_cft\n"
            + "r = prep_axloc(%d)\n" % c
            + "plt.figure()\n"
            + "if r is None:\n"
            + "    print('no diluted local axial')\n"
            + "else:\n"
            + "    m,ere,eim,d = r\n"
            + "    v = m%s; e = %s; dv = d%s\n" % (suf, err, suf)
            + "    pos = v[dt]>0; neg = v[dt]<0; dpos = dv[dt]>0; dneg = dv[dt]<0\n"
            + "    plt.errorbar(t[pos], v[dt][pos], yerr=e[dt][pos], fmt='o', ms=3, capsize=2, label='diluted $+$')\n"
            + "    plt.errorbar(t[neg], -v[dt][neg], yerr=e[dt][neg], fmt='v', ms=3, capsize=2, mfc='none', label='diluted $-$')\n"
            + "    plt.plot(t[dpos], dv[dt][dpos], 's--', ms=3, label='determ $+$')\n"
            + "    plt.plot(t[dneg], -dv[dt][dneg], 'x:', ms=4, label='determ $-$')\n"
            + "    plt.yscale('log'); plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$|\\mathrm{%s}|$ local axial $s_%d$')\n" % (P, c)
            + "    plt.title('axial local $s_%d$ %s (signed log), $m='+MTIT+'$'); plt.legend(fontsize=8)\n" % (c, P)
            + "    plt.savefig(FIGDIR+'/ax_local_s%d_%s_'+MLAB+'.pdf', bbox_inches='tight')\n" % (c, part))

# ---------------------------------------------------------------- non-split cells (ratios/summary/variance/condensate)
VEC_HEADER = r'''## Vector (conserved-current $J$)
Exact-K temporal $G_t$ (`tp`) and spatial $G_s$ (`sp`), plus the local $\sigma_a$ current.  Diluted (jk,
Nt/2-shift averaged) vs the determ ground truth.  LOG axis; the $+$ and $-$ sign branches of each Re/Im part
are in separate cells (the branch with data shows the decay).  $0<t<N_t/2$.'''

VEC_RATIO = r'''# VECTOR exact ratio Gs/Gt (Re; normalization-free).  diluted (jk) vs determ.  LINEAR.  0<t<Nt/2.
dt = np.arange(1, Nt_cft//2); t = dt*at_cft; ref = 5
sp = half_shift_avg(dil_hits(DIL,'sp/Vpp')); tp = half_shift_avg(dil_hits(DIL,'tp/Vpp'))
m,e = jk_ratio_re(sp, tp)
plt.figure()
plt.errorbar(t, m[dt], yerr=e[dt], fmt='o', ms=3, capsize=2, label='diluted (jk)')
ds = half_shift_avg(mirror(det_raw(EXS,'sp/Vpp'))).real; dtt = half_shift_avg(mirror(det_raw(EXS,'tp/Vpp'))).real
plt.plot(t, (ds/dtt)[dt], 's--', ms=3, label='determ exact-K')
plt.axhline(-1.0,color='k',lw=1,label='CFT $-1$'); plt.ylim(-3,0.5)
plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$G_s/G_t$')
plt.title('vector exact ratio $G_s/G_t$, $m='+MTIT+'$'); plt.legend(fontsize=8)
plt.savefig(FIGDIR+'/vec_exact_ratio_'+MLAB+'.pdf', bbox_inches='tight')
print('exact ratio @dt=%d: % .3f +- %.3f'%(ref,m[ref],e[ref]))
'''

VEC_LOC_RATIO = r'''# VECTOR local ratio (s1+s2)/2 / s3 -> -1 (Re).  diluted (jk) vs determ.  LINEAR.  0<t<Nt/2.
dt = np.arange(1, Nt_cft//2); t = dt*at_cft; ref = 5
num = 0.5*(half_shift_avg(dil_hits(DIL,'s1/Vpp'))+half_shift_avg(dil_hits(DIL,'s2/Vpp')))
den = half_shift_avg(dil_hits(DIL,'s3/Vpp'))
m,e = jk_ratio_re(num, den)
plt.figure()
plt.errorbar(t, m[dt], yerr=e[dt], fmt='o', ms=3, capsize=2, label='diluted (jk)')
dn = 0.5*(half_shift_avg(mirror(det_raw(LOC,'s1/Vpp')))+half_shift_avg(mirror(det_raw(LOC,'s2/Vpp')))).real
dd = half_shift_avg(mirror(det_raw(LOC,'s3/Vpp'))).real
plt.plot(t, (dn/dd)[dt], 's--', ms=3, label='determ local')
plt.axhline(-1.0,color='k',lw=1,label='CFT $-1$'); plt.ylim(-2,0.5)
plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$(s_1{+}s_2)/(2 s_3)$')
plt.title('vector local ratio, $m='+MTIT+'$'); plt.legend(fontsize=8)
plt.savefig(FIGDIR+'/vec_local_ratio_'+MLAB+'.pdf', bbox_inches='tight')
print('local ratio @dt=%d: % .3f +- %.3f'%(ref,m[ref],e[ref]))
'''

VEC_LOC_S1S2 = r'''# VECTOR local |s1| vs |s2| (free: spin-diagonal => s1==s2).  0<t<Nt/2, log.
dt = np.arange(1, Nt_cft//2); t = dt*at_cft
m1,_,_ = jk_cplx(half_shift_avg(dil_hits(DIL,'s1/Vpp'))); m2,_,_ = jk_cplx(half_shift_avg(dil_hits(DIL,'s2/Vpp')))
plt.figure()
plt.plot(t, np.abs(m1[dt]),'o',ms=3,label='$|s_1|$')
plt.plot(t, np.abs(m2[dt]),'x',ms=4,label='$|s_2|$')
plt.yscale('log'); plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$|s_a|$')
plt.title('vector local $s_1$ vs $s_2$, $m='+MTIT+'$'); plt.legend(fontsize=8)
plt.savefig(FIGDIR+'/vec_local_s1s2_'+MLAB+'.pdf', bbox_inches='tight')
print('max|s1-s2|/|s1| over dt:', np.max(np.abs(m1[dt]-m2[dt])/np.abs(m1[dt])))
'''

VEC_VARIANCE = r'''# VARIANCE: diluted vs volume error bars on -Re G_t (Nt/2-shift avg, log).  Volume run exists only for m0;
# this cell auto-skips for mF/mP.  Matched cost ~ 1 diluted hit = 4 patterns ~ 4 volume hits.
dt = np.arange(1, Nt_cft//2); t = dt*at_cft
plt.figure()
m,ere,eim,d = prep_v('tp','Vpp')
plt.errorbar(t, -m.real[dt], yerr=ere[dt], fmt='o', ms=3, capsize=2, label='diluted (%d hits x4 patt)'%len(dil_files(DIL)))
vfs = [f for f in sorted(glob.glob(ESN+VOL+'/corr.0.h*.h5')) if ('complete' in h5py.File(f,'r'))]
if vfs:
    vh=[]
    for f in vfs:
        with h5py.File(f,'r') as h:
            bs=[k for k in h['h0'] if k.startswith('t0_')]; acc=None
            for b in bs:
                v=h['h0/'+b+'/tp/Vpp/real'][()]+1j*h['h0/'+b+'/tp/Vpp/imag'][()]
                acc=v if acc is None else acc+v
            vh.append(acc/len(bs))
    mv,ev,_ = jk_cplx(half_shift_avg(np.array(vh)))
    plt.errorbar(t+0.01, -mv.real[dt], yerr=ev[dt], fmt='s', ms=3, capsize=2, label='volume (%d hits)'%len(vh))
else:
    print('volume run not found (only m0 has it) -- diluted shown alone')
plt.yscale('log'); plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$-\mathrm{Re}\,G_t$')
plt.title('variance: diluted vs volume, $m='+MTIT+'$'); plt.legend(fontsize=8)
plt.savefig(FIGDIR+'/vec_variance_'+MLAB+'.pdf', bbox_inches='tight')
'''

VEC_SUMMARY = r'''# VECTOR summary @dt=5 (Nt/2-shift avg, Re).  diluted vs determ.
ref = 5
mt,et,_ = jk_cplx(half_shift_avg(dil_hits(DIL,'tp/Vpp'))); ms,es,_ = jk_cplx(half_shift_avg(dil_hits(DIL,'sp/Vpp')))
mr,er = jk_ratio_re(half_shift_avg(dil_hits(DIL,'sp/Vpp')), half_shift_avg(dil_hits(DIL,'tp/Vpp')))
print('EXACT  @dt=%d:  G_t=% .3e +- %.1e   G_s=% .3e +- %.1e   G_s/G_t=% .3f +- %.3f'%(ref,mt.real[ref],et[ref],ms.real[ref],es[ref],mr[ref],er[ref]))
dgt = half_shift_avg(mirror(EXS_SCALE*det_raw(EXS,'tp/Vpp'))).real; dgs = half_shift_avg(mirror(EXS_SCALE*det_raw(EXS,'sp/Vpp'))).real
print('DETERM @dt=%d:  G_t=% .3e   G_s=% .3e   G_s/G_t=% .3f'%(ref,dgt[ref],dgs[ref],dgs[ref]/dgt[ref]))
'''

AX_HEADER = r'''## Axial current ($C_{A+-}$)
Exact-K axial `tp`/`sp` (`h0/axial/*/Apm`) and the local axial $s_3$.  Diluted (jk, $\times c=1/4\pi$) vs
determ exact1.  LOG axis; $+$ and $-$ sign branches of each Re/Im part in separate cells.  $0<t<N_t/2$.'''

AX_LOC_RATIO = r'''# AXIAL local ratio (s1+s2)/2 / s3 (Re).  diluted (jk) vs determ.  LINEAR.  0<t<Nt/2.
dt = np.arange(1, Nt_cft//2); t = dt*at_cft
A1,A2,A3 = dil_ax_loc(1), dil_ax_loc(2), dil_ax_loc(3)
plt.figure()
if A1 is None:
    print('no diluted local axial')
else:
    rn,_,_ = jk_cplx(0.5*(half_shift_avg(A1)+half_shift_avg(A2))); rd,_,_ = jk_cplx(half_shift_avg(A3))
    plt.plot(t, (rn.real/rd.real)[dt], 'o', ms=3, label='diluted local')
    d1 = half_shift_avg(mirror(det_raw(LOCAX,'s1/Apm'))); d2 = half_shift_avg(mirror(det_raw(LOCAX,'s2/Apm'))); d3 = half_shift_avg(mirror(det_raw(LOCAX,'s3/Apm')))
    plt.plot(t, ((0.5*(d1+d2)).real/d3.real)[dt], 's--', ms=3, label='determ local axial')
    plt.axhline(-1.0,color='k',lw=1,label='$-1$'); plt.ylim(-3,1)
    plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$(s_1{+}s_2)/(2 s_3)$ local axial')
    plt.title('axial local ratio, $m='+MTIT+'$'); plt.legend(fontsize=8)
    plt.savefig(FIGDIR+'/ax_local_ratio_'+MLAB+'.pdf', bbox_inches='tight')
'''

AX_SUMMARY = r'''# AXIAL summary: diluted (x c) vs determ, Re and Im, both channels, a few dt.
for proj in ['tp','sp']:
    m,ere,eim,d,c = prep_ax(proj)
    print('=== axial %s (c=%.4f) ===' % (proj, c))
    for dt_ in [1,3,5,8,12]:
        print(' dt%2d  Re dil % .3e+-%.1e det % .3e | Im dil % .3e+-%.1e det % .3e' % (
            dt_, (c*m).real[dt_], c*ere[dt_], d.real[dt_], (c*m).imag[dt_], c*eim[dt_], d.imag[dt_]))
'''

COND_HEADER = r'''## Condensates ($\sigma_{PS}$ Eq. 1.23, $\sigma_{FS}$ Eq. 1.55)
Spacetime-averaged ($/V_{st}$) bilinears from `h0/condensate/*` (jk over hits) vs the DENSE exact reference
(`condensate_deter_L1`).  $\sigma_{PS}=$`etadag_xi`$+$`xidag_eta`, $\sigma_{FS}=$`etadag_xi`$-$`xidag_1mDdag_eta`.'''

COND_LOAD = r'''# ----- condensate: dense ref + stochastic (jk) + table -----
DENSE = ESN + 'condensate_deter_L1/cond.h5'
STOCH = ESN + DIL + '/'
CKEYS = ['etadag_xi', 'xidag_eta', 'xidag_1mDdag_eta']
def cond_stoch(key):
    out=[]
    for fn in sorted(glob.glob(STOCH+'corr.0.h*.h5')):
        with h5py.File(fn,'r') as h:
            if 'complete' not in h: continue
            kk='h0/condensate/'+key
            out.append(h[kk+'/real'][()][0] + 1j*h[kk+'/imag'][()][0])
    return np.array(out)
def jk_scalar(a):
    H=len(a); s=a.sum()
    samp=np.array([(s-a[i])/(H-1) for i in range(H)])
    est=samp.mean(); vre=(H-1)*np.mean((samp.real-est.real)**2); vim=(H-1)*np.mean((samp.imag-est.imag)**2)
    return est, np.sqrt(vre), np.sqrt(vim)
with h5py.File(DENSE,'r') as f:
    Vst = f['Vst'][()][0]
    dense = {k: f['condensate/'+k+'/real'][()][0] + 1j*f['condensate/'+k+'/imag'][()][0] for k in CKEYS}
dense['sigma_PS'] = dense['etadag_xi'] + dense['xidag_eta']
dense['sigma_FS'] = dense['etadag_xi'] - dense['xidag_1mDdag_eta']
raw = {k: cond_stoch(k) for k in CKEYS}
nh = len(raw['etadag_xi']); assert nh > 0, 'no condensate hits in %s' % STOCH
raw['sigma_PS'] = raw['etadag_xi'] + raw['xidag_eta']
raw['sigma_FS'] = raw['etadag_xi'] - raw['xidag_1mDdag_eta']
print('condensate hits: %d   Vst=%.4f' % (nh, Vst))
print('%-18s %26s %26s' % ('quantity','stoch /Vst (re,im)','dense /Vst (re,im)'))
for k in CKEYS + ['sigma_PS','sigma_FS']:
    e,er,ei = jk_scalar(raw[k]); e,er,ei = e/Vst,er/Vst,ei/Vst; dd = dense[k]/Vst
    print('%-18s % .4e%+.4ej(+-%.0e,%.0e) % .4e%+.4ej' % (k,e.real,e.imag,er,ei,dd.real,dd.imag))
'''

def cond_plot(sig, eq):
    tex = sig.replace('_', '_{') + '}'   # sigma_PS -> sigma_{PS} (both letters subscripted)
    # exact momentum-independent (delta_xy) contact per channel/run -- condensate_contact_massive_claude.md Sec 10:
    if sig == 'sigma_PS':
        csub_expr = "(1.0+1.0/(1.0+m_cur)) if MASS=='mP' else 2.0"   # PS: 2 (m0,mF); 1+1/(1+m_P) (mP)
    else:  # sigma_FS: 2 (m0); 2-m_F (mF, clean dagger bwd); 2-(m_P/(1+m_P))^2 (mP, rescaled bwd)
        csub_expr = "2.0 if MASS=='massless' else (2.0-m_cur if MASS=='mF' else 2.0-(m_cur/(1.0+m_cur))**2)"
    return ("# CONDENSATE %s, CONTACT-SUBTRACTED (condensate_contact_massive_claude.md Sec 10: exact m.i. contact).\n"
            "# contact of tr[A X] in /Vst units = 2*[[X]]; forward leg 1/2, backward 1/2 (dagger) or 1/(2(1+m_P)).\n"
            "# Dynamical remainder ~0 except sigma_PS(m_F) / sigma_FS(m_P).\n"
            "e,er,ei = jk_scalar(raw['%s'])\n"
            "m_cur = {'massless':0.0,'mF':0.1+0j,'mP':0.1j}[MASS]; csub = %s\n"
            "es = e/Vst + csub; dd = dense['%s']/Vst + csub\n"
            "plt.figure()\n"
            "plt.errorbar([0,1],[es.real,es.imag],yerr=[er/Vst,ei/Vst],fmt='o',capsize=4,label='stoch (jk)')\n"
            "plt.plot([0,1],[dd.real,dd.imag],'s',ms=9,mfc='none',label='dense exact')\n"
            "plt.axhline(0,color='gray',lw=0.5); plt.xticks([0,1],['Re','Im'])\n"
            "plt.ylabel(r'$\\langle\\%s\\rangle/V_{st}-$contact'); plt.title(r'$\\%s$ (%s) contact-subtracted, $m='+MTIT+'$')\n"
            "plt.legend(); plt.grid(alpha=0.3)\n"
            "plt.savefig(FIGDIR+'/cond_%s_sub_'+MLAB+'.pdf', bbox_inches='tight')\n"
            % (sig, csub_expr, sig, tex, tex, eq, sig))

# ---------------------------------------------------------------- DISC + physical-sum cells
DISC_HEADER = r'''## Disconnected current ($C_\mathrm{disc}$) and physical sum
Disc post-processing (mirrors `jj_disc_postproc_claude.cc`; the C++ script reads the OLD multi-hit-per-file
layout, so for the dilute per-hit files it is done inline here): the disc one-point trace $J(t)$ is RAW
(no $1/4\pi$) per hit, and
$$C^P_\mathrm{disc}(\Delta t)=\tfrac{1}{4\pi}\,\tfrac{1}{N_t}\sum_{t_0}\bar J(t_0)\,\bar J(t_0+\Delta t),\quad
\bar J=\text{hit-average of }J.$$
For the FREE case ($U=1$, single config) the disc correlator is $\approx 0$ (noise floor $\sim1/N_\mathrm{hits}$).
The physical correlator is then $-C_\mathrm{conn}+C_\mathrm{disc}$ (PDF Eq. 3.39, vector only; conn $=V_{++}$
which carries its own $1/4\pi$).  Signed log, both branches, $0<t<N_t/2$.'''

DISC_LOAD = r'''# ----- disc post-processing (mirrors jj_disc_postproc_claude.cc) -----
INV4PI = 1.0/(4.0*np.pi)
def disc_J(proj, til=False):             # per-hit RAW one-point trace -> (nhits, Nt) complex
    key = 'h0/disc/%s/%s' % (proj, 'Jtil' if til else 'J')
    out=[]
    for f in dil_files(DIL):
        with h5py.File(f,'r') as h:
            out.append(h[key+'/real'][()] + 1j*h[key+'/imag'][()])
    return np.array(out)
def two_point(A, B):                     # (1/Nt) sum_t0 A[t0] B[(t0+dt) mod Nt] -> (Nt,) complex
    Nt=len(A); return np.array([np.mean(A*np.roll(B,-dt)) for dt in range(Nt)])
def prep_disc(proj, pm='Vpp'):           # disc two-point (jk over hits), Nt/2-shift'd -> (m, ere, eim)
    til = (pm=='Vmm' and parity); J = disc_J(proj, til=til); H=len(J); S=J.sum(0)
    samp=[]
    for i in range(H):
        jb=(S-J[i])/(H-1); c=INV4PI*two_point(jb,jb)
        if pm=='Vmm' and not parity: c=np.conj(c)
        samp.append(half_shift_avg(c))
    samp=np.array(samp); est=samp.mean(0)
    ere=np.sqrt((H-1)*np.mean((samp.real-est.real)**2,0)); eim=np.sqrt((H-1)*np.mean((samp.imag-est.imag)**2,0))
    return est, ere, eim
def prep_phys(proj):                     # physical = -C_conn + C_disc (Eq 3.39, vector V++), Nt/2-shift'd.
    # SINGLE delete-1 jackknife over hits: each pseudo-value recomputes BOTH the conn mean and the
    # (quadratic) disc two-point on the SAME leave-one-out hit set, so the error carries the disc variance
    # AND the conn-disc covariance (disc is NOT added as a fixed constant -- that undercounts the error).
    conn = half_shift_avg(dil_hits(DIL, proj+'/Vpp')); J = disc_J(proj)
    H = len(conn); Sc = conn.sum(0); SJ = J.sum(0)
    samp = []
    for i in range(H):
        Jbar_i = (SJ - J[i])/(H-1)
        samp.append(-(Sc - conn[i])/(H-1) + half_shift_avg(INV4PI*two_point(Jbar_i, Jbar_i)))
    samp = np.array(samp); est = samp.mean(0)
    ere = np.sqrt((H-1)*np.mean((samp.real-est.real)**2,0))
    eim = np.sqrt((H-1)*np.mean((samp.imag-est.imag)**2,0))
    d = -half_shift_avg(mirror(EXS_SCALE*det_raw(EXS, proj+'/Vpp')))   # -conn determ (+ disc determ = 0 free)
    return est, ere, eim, d
print('disc one-point hits:', len(dil_files(DIL)))
'''

def disc_cell(proj, part):               # disc-only signed-log (diluted only; free -> ~0)
    P   = 'Re' if part=='re' else 'Im'
    suf = '.real' if part=='re' else '.imag'
    err = 'ere' if part=='re' else 'eim'
    G   = 't' if proj=='tp' else 's'
    return ("# DISC two-point %s V++ %s on signed-log (free -> ~0 noise floor).  0<t<Nt/2.\n" % (proj, P)
            + "dt = np.arange(1, Nt_cft//2); t = dt*at_cft\n"
            + "m,ere,eim = prep_disc('%s')\n" % proj
            + "v = m%s; e = %s\n" % (suf, err)
            + "pos = v[dt] > 0; neg = v[dt] < 0\n"
            + "plt.figure()\n"
            + "plt.errorbar(t[pos], v[dt][pos], yerr=e[dt][pos], fmt='o', ms=3, capsize=2, label='disc $+$')\n"
            + "plt.errorbar(t[neg], -v[dt][neg], yerr=e[dt][neg], fmt='v', ms=3, capsize=2, mfc='none', label='disc $-$')\n"
            + "plt.yscale('log'); plt.xlabel(r'$t=a_t n_t$'); plt.ylabel(r'$|\\mathrm{%s}\\,C^{disc}_{%s}|$')\n" % (P, G)
            + "plt.title('disc %s %s (free $\\\\to 0$), $m='+MTIT+'$'); plt.legend(fontsize=8)\n" % (proj, P)
            + "plt.savefig(FIGDIR+'/disc_%s_%s_'+MLAB+'.pdf', bbox_inches='tight')\n" % (proj, part))

SUM_HEADER = r'''### Physical sum $-C_\mathrm{conn}+C_\mathrm{disc}$ (Eq. 3.39, vector)
Conn $=V_{++}$ exact-K (diluted); disc $\approx0$ for free, so the physical correlator $\approx -C_\mathrm{conn}$
(sign flips relative to conn: $G_t$ becomes positive, $G_s$ negative).  Determ overlay $=-C_\mathrm{conn}^\mathrm{determ}$.'''

TITLE = {
    'massless': "# Validation -- massless ($m=0$)\nVector + axial + condensate vs determ/dense ground truth (free field, $L=1$).  Nt/2-shift averaged; log-scale correlators with +/- sign branches in separate cells; $0<t<N_t/2$.",
    'mF':       "# Validation -- $m_F=0.1$ (flavor-breaking, real mass)\nVector + axial + condensate vs determ/dense ground truth (free field, $L=1$).  Nt/2-shift averaged; log-scale correlators with +/- sign branches in separate cells; $0<t<N_t/2$.",
    'mP':       "# Validation -- $m_P=0.1i$ (parity-breaking, imaginary mass)\nVector + axial + condensate vs determ/dense ground truth (free field, $L=1$).  $V_{++}$ and $V_{--}$ are independent here.  Nt/2-shift averaged; log-scale correlators with +/- sign branches in separate cells; $0<t<N_t/2$.",
}

# prep statements (set m,ere,eim,d[,c]) shared by the +/- cell pairs
def prep_stmt(call):
    return "dt = np.arange(1, Nt_cft//2); t = dt*at_cft\n" + call + "\n"

def build_cells(mass):
    C = []
    def md(s):   C.append(('markdown', s))
    def code(s): C.append(('code', s))
    md(TITLE[mass]); code(CONFIG.replace('__MASS__', mass)); code(HELPERS)
    # ---- vector ----
    md(VEC_HEADER)
    for proj in ['tp','sp']:
        G = 'G_t' if proj=='tp' else 'G_s'
        prep = prep_stmt("m,ere,eim,d = prep_v('%s','Vpp')" % proj)
        for part in ['re','im']:
            code(cell(prep, part, '%s^{++}'%G, 'vector exact $%s$ $V_{++}$'%G, 'vec_exact_%s_pp'%proj))
    code(VEC_RATIO)
    if mass == 'mP':
        for proj in ['tp','sp']:
            G = 'G_t' if proj=='tp' else 'G_s'
            prep = prep_stmt("m,ere,eim,d = prep_v('%s','Vmm')" % proj)
            for part in ['re','im']:
                code(cell(prep, part, '%s^{--}'%G, 'vector exact $%s$ $V_{--}$'%G, 'vec_exact_%s_mm'%proj))
    prep = prep_stmt("m,ere,eim,d = prep_loc(3)")
    for part in ['re','im']:
        code(cell(prep, part, 's_3', 'vector local $s_3$', 'vec_local_s3'))
    prep = prep_stmt("m,ere,eim,d = prep_locsp()")
    for part in ['re','im']:
        code(cell(prep, part, '(s_1{+}s_2)/2', 'vector local $G_s$', 'vec_local_sp'))
    code(VEC_LOC_RATIO); code(VEC_LOC_S1S2); code(VEC_VARIANCE); code(VEC_SUMMARY)
    # ---- axial ----
    md(AX_HEADER)
    for proj in ['tp','sp']:
        prep = prep_stmt("m,ere,eim,d,c = prep_ax('%s')" % proj)
        for part in ['re','im']:
            code(cell(prep, part, 'C_{A+-}^{%s}'%proj, 'axial %s'%proj, 'ax_exact_%s'%proj, cfac=True))
    for part in ['re','im']:
        code(axloc_cell(3, part))
    code(AX_LOC_RATIO); code(AX_SUMMARY)
    # ---- condensate ----
    md(COND_HEADER); code(COND_LOAD); code(cond_plot('sigma_PS','1.23')); code(cond_plot('sigma_FS','1.55'))
    # ---- disc + physical sum ----
    md(DISC_HEADER); code(DISC_LOAD)
    for proj in ['tp','sp']:
        for part in ['re','im']:
            code(disc_cell(proj, part))
    md(SUM_HEADER)
    for proj in ['tp','sp']:
        G = 'G_t' if proj=='tp' else 'G_s'
        prep = prep_stmt("m,ere,eim,d = prep_phys('%s')" % proj)
        for part in ['re','im']:
            code(cell(prep, part, '%s'%G, 'physical $-C_c{+}C_d$ $%s$'%G, 'phys_%s'%proj))
    return C

def to_nb(cells):
    out=[]
    for ctype, src in cells:
        cell={'cell_type':ctype,'metadata':{},'source':src,'id':uuid.uuid4().hex[:12]}
        if ctype=='code': cell['execution_count']=None; cell['outputs']=[]
        out.append(cell)
    return {'cells':out,
            'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},
                        'language_info':{'name':'python','version':'3'}},
            'nbformat':4,'nbformat_minor':5}

for mass, fn in [('massless','jj_validate_m0_claude.ipynb'),
                 ('mF','jj_validate_mF_claude.ipynb'),
                 ('mP','jj_validate_mP_claude.ipynb')]:
    nb = to_nb(build_cells(mass))
    with open(fn,'w') as f:
        json.dump(nb, f, indent=1)
    print('wrote %s  (%d cells)' % (fn, len(nb['cells'])))
