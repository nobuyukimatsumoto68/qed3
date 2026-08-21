# effmass_axial_sp_expfit_claude.py  (sp = SO(3) spatial scalar s1+s2; cf tp = s3)
# Const+exp effmass fit  m_eff(dt) = m0 + A exp(-c1 dt)  of the AXIAL-current SP (spatial, s1+s2) l=1 (m-averaged,
# conn) cosh effmass.  m0 = plateau mass.  Fit windows (dt = timeslice sep, NM 2026-08-18):
#   L1, L2: dt in [10,25]   ;   L3, L4: dt in [5,15]
# Config-jackknife (fit each delete-1 sample -> m0^b).  kmin=20.  Nf 2/4/6 overlaid; L1..L4 panels.
# The PNG overlays the fitted m0+A exp(-c1 dt) curve on the effmass points.
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
# axial tp const+exp fit windows (dt = timeslice sep, inclusive), NM 2026-08-18.
# BASE per L: L1/L2 [10,25], L3/L4 [5,15].  PER-(ell,L) OVERRIDES below.
FITDT_BASE = {1: (10, 25), 2: (10, 25), 3: (5, 15), 4: (5, 15)}
# overrides keyed (ell, L): ell=1 at L3,L4 -> [5,20] (NM 2026-08-18)
FITDT_OVR = {(1, 3): (5, 20), (1, 4): (5, 20)}


def win(ell, L):
    lo, hi = FITDT_OVR.get((ell, L), FITDT_BASE[L])
    return np.arange(lo, hi + 1)


dtp = np.arange(1, Nt // 2)
ELLS = [0, 1]                          # sp: ell=0 and ell=1 (NM 2026-08-18: ell=2,3 not needed for sp)
ELL_L = {0: [1, 2, 3, 4], 1: [1, 2, 3, 4], 2: [1, 2, 3, 4], 3: [1, 2]}


def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB[L]))


def conn_files(nf, L, g):
    fs = sorted(glob.glob(esn(nf, L, g) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if int(re.search(r"corr\.(\d+)\.", f).group(1)) >= KMIN]


def load(fs, key):
    o = []
    for fn in fs:
        with h5py.File(fn, 'r') as f:
            o.append(f[key + '/real'][()] + 1j * f[key + '/imag'][()])
    return np.array(o)


def C_axial(fs, l):
    # sp = SO(3) spatial scalar = s1 + s2 (m-summed), per jj_ylm_prod_analysis_redo_claude.ipynb.
    ssum = sum(load(fs, "h0/ylm_axial/s1/l%d/m%d/Vpp" % (l, m)) + load(fs, "h0/ylm_axial/s2/l%d/m%d/Vpp" % (l, m))
               for m in range(-l, l + 1)) / (2 * l + 1)
    return 2.0 * ssum.real


def meff_acosh(C):
    m = np.full(Nt, np.nan)
    for t in range(1, Nt - 1):
        d = 2.0 * C[t]
        if d == 0:
            continue
        r = (C[t - 1] + C[t + 1]) / d
        if r > 1.0:
            m[t] = np.arccosh(r)
    return m


def effmass_jk(samp):
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def model(dt, m0, A, c1):
    return m0 + A * np.exp(-c1 * dt)


def m0fit(me, mean, err, D):
    x = D.astype(float)
    if np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
        return None, None, None
    p0 = [mean[D[-1]], mean[D[0]] - mean[D[-1]], 0.3]
    bnd = ([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf])
    try:
        popt, _ = curve_fit(model, x, mean[D], sigma=err[D], p0=p0, bounds=bnd, maxfev=40000)
    except (RuntimeError, ValueError):
        return None, None, None
    H = me.shape[0]
    m0b = np.full(H, np.nan)
    for b in range(H):
        yb = me[b, D]
        if np.any(~np.isfinite(yb)):
            continue
        try:
            pb, _ = curve_fit(model, x, yb, sigma=err[D], p0=popt, bounds=bnd, maxfev=40000)
            m0b[b] = pb[0]
        except (RuntimeError, ValueError):
            continue
    good = np.isfinite(m0b)
    if good.sum() < H // 2:
        return None, None, None
    mb = m0b[good]
    n = mb.size
    M = mb.mean()
    sig = math.sqrt((n - 1.0) / n * np.sum((mb - M) ** 2))
    return M, sig, popt


# ---------- axial tp const+exp fit + effmass plot, ONE png per ell (fit curve overlaid) ----------
gls = {0.5: "-", 1.0: "--", 1.5: ":", 2.0: "--", 3.0: ":", 4.5: "-.", 4.0: "--", 6.0: "-."}
irrep = {0: "A", 1: "T1", 2: "H", 3: "T2+G"}
ymax = {0: 4, 1: 4, 2: 5, 3: 6}   # higher ell = heavier -> more head-room

for LVAL in ELLS:
    rows = []   # (L, g, nf, m0, sig, ncfg)
    Ls = ELL_L[LVAL]
    fig, axs = plt.subplots(1, len(Ls), figsize=(6 * len(Ls), 5))
    for ax, L in zip(np.atleast_1d(axs).ravel(), Ls):
        D = win(LVAL, L)
        xc = np.linspace(D[0], D[-1], 60)
        for g in GS[L]:
            for nf in NFS:
                fs = conn_files(nf, L, g)
                if len(fs) < 4:
                    continue
                me, mean, err = effmass_jk(C_axial(fs, LVAL).astype(float))
                M, sig, popt = m0fit(me, mean, err, D)
                ax.errorbar(dtp, mean[dtp], yerr=err[dtp], fmt='o', ms=2.5, capsize=1.5,
                            color=nfcol[nf], ls="none", alpha=0.6, label="Nf%d g%.1f" % (nf, g))
                if M is not None:
                    ax.plot(xc, model(xc, *popt), color=nfcol[nf], lw=1.6, ls=gls[g])   # FIT CURVE
                    ax.axhline(M, color=nfcol[nf], ls=gls[g], lw=0.7, alpha=0.5)          # m0 plateau
                    rows.append((L, g, nf, M, sig, len(fs)))
        ax.axvspan(D[0], D[-1], color="gray", alpha=0.1)
        ax.set_title(r"axial sp $\ell=%d$ (%s, m-avg) L%d  fit dt[%d,%d]"
                     % (LVAL, irrep[LVAL], L, D[0], D[-1]))
        ax.set_xlim(0, 30)
        ax.set_ylim(0, ymax[LVAL])
        ax.set_xlabel("dt")
        ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6, ncol=3)
    fig.suptitle(r"AXIAL sp $\ell=%d$ (%s): fit $m_0+A\,e^{-c_1 dt}$  (m0 = plateau; curve overlaid)"
                 % (LVAL, irrep[LVAL]))
    fig.tight_layout()
    fig.savefig("figs/effmass_axial_sp_l%d_expfit_claude.png" % LVAL, dpi=150)
    print("# wrote effmass_axial_tp_l%d_expfit_claude.png (%d fits)" % (LVAL, len(rows)))

    out = ["# Axial sp l=%d (%s) const+exp fit masses (m-avg conn, m0+A exp(-c1 dt), kmin=20)"
           % (LVAL, irrep[LVAL]), "",
           "a_t m0 from m_eff(dt)=m0+A exp(-c1 dt), config-jackknife.",
           "windows (dt): L1/L2 [10,25], L3/L4 [5,15].", "",
           "| L | gsq | Nf | a_t m0 | ncfg |", "|---|-----|----|--------|------|"]
    for L, g, nf, M, sig, n in rows:
        out.append("| %d | %.1f | %d | %.4f(%.4f) | %d |" % (L, g, nf, M, sig, n))
    open("axial_sp_l%d_expfit_masses_claude.md" % LVAL, "w").write("\n".join(out) + "\n")
