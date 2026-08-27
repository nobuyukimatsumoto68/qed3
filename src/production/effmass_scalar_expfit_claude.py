# effmass_scalar_expfit_claude.py
# Const+exp effmass fit  m_eff(dt) = m0 + A exp(-c1 dt)  of the SCALAR conn correlators, BOTH channels:
#   PS = 2*Vpp                    (parity-symmetric)
#   FS = Vpp + Vmm^FS             (furnished; scalar_fs Vmm)
# m-averaged, config-jackknife, cosh effmass, kmin=20.  Fit curve overlaid on the effmass points.
# Ells / windows (dt, inclusive), NM 2026-08-18:
#   ell=1: L1,L2 [8,24] ; L3,L4 [6,16]
#   ell=2: L1,L2 [12,22]   (ell=2 only for L1,L2; L3,L4 use ell=1 only)
# One png per (channel, ell): effmass_scalar_<PS|FS>_l<ell>_expfit_claude.png (L panels side by side).
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
CHANS = ["PS", "FS"]
# which L per ell:
ELL_L = {1: [1, 2, 3, 4], 2: [1, 2]}
# per-(ell,L) fit windows (dt, inclusive)
WIN = {(1, 1): (8, 24), (1, 2): (8, 24), (1, 3): (6, 16), (1, 4): (6, 16),
       (2, 1): (12, 22), (2, 2): (10, 16)}
dtp = np.arange(1, Nt // 2)


def win(ell, L):
    lo, hi = WIN[(ell, L)]
    return np.arange(lo, hi + 1)


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


def C_scalar(fs, l, chan):
    vpp = sum(load(fs, "h0/scalar/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)
    if chan == "PS":
        return (2.0 * vpp).real
    vfs = sum(load(fs, "h0/scalar_fs/l%d/m%d/Vmm" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)
    return (vpp + vfs).real


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


gls = {0.5: "-", 1.0: "--", 1.5: ":", 2.0: "--", 3.0: ":", 4.5: "-.", 4.0: "--", 6.0: "-."}
ymax = {1: 4, 2: 5}

for chan in CHANS:
    for ell in (1, 2):
        Ls = ELL_L[ell]
        rows = []
        fig, axs = plt.subplots(1, len(Ls), figsize=(6.5 * len(Ls), 5))
        for ax, L in zip(np.atleast_1d(axs).ravel(), Ls):
            D = win(ell, L)
            xc = np.linspace(D[0], D[-1], 60)
            for g in GS[L]:
                for nf in NFS:
                    fs = conn_files(nf, L, g)
                    if len(fs) < 4:
                        continue
                    me, mean, err = effmass_jk(C_scalar(fs, ell, chan).astype(float))
                    M, sig, popt = m0fit(me, mean, err, D)
                    ax.errorbar(dtp, mean[dtp], yerr=err[dtp], fmt='o', ms=2.5, capsize=1.5,
                                color=nfcol[nf], ls="none", alpha=0.6, label="Nf%d g%.1f" % (nf, g))
                    if M is not None:
                        ax.plot(xc, model(xc, *popt), color=nfcol[nf], lw=1.6, ls=gls[g])
                        ax.axhline(M, color=nfcol[nf], ls=gls[g], lw=0.7, alpha=0.5)
                        rows.append((L, g, nf, M, sig, len(fs)))
            ax.axvspan(D[0], D[-1], color="gray", alpha=0.1)
            ax.set_title(r"scalar %s $\ell=%d$ L%d  fit dt[%d,%d]" % (chan, ell, L, D[0], D[-1]))
            ax.set_xlim(0, 30)
            ax.set_ylim(0, ymax[ell])
            ax.set_xlabel("dt")
            ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=6, ncol=3)
        fig.suptitle(r"SCALAR %s $\ell=%d$: fit $m_0+A\,e^{-c_1 dt}$  (m0 = plateau; curve overlaid)"
                     % (chan, ell))
        fig.tight_layout()
        fig.savefig("figs/effmass_scalar_%s_l%d_expfit_claude.png" % (chan, ell), dpi=150)
        print("# wrote effmass_scalar_%s_l%d_expfit_claude.png (%d fits)" % (chan, ell, len(rows)))

        out = ["# Scalar %s l=%d const+exp fit masses (m-avg conn, m0+A exp(-c1 dt), kmin=20)"
               % (chan, ell), "",
               "a_t m0 from m_eff(dt)=m0+A exp(-c1 dt), config-jackknife.", "",
               "| L | gsq | Nf | a_t m0 | ncfg |", "|---|-----|----|--------|------|"]
        for L, g, nf, M, sig, n in rows:
            out.append("| %d | %.1f | %d | %.4f(%.4f) | %d |" % (L, g, nf, M, sig, n))
        open("scalar_%s_l%d_expfit_masses_claude.md" % (chan, ell), "w").write("\n".join(out) + "\n")
