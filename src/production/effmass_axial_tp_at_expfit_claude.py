# effmass_axial_tp_at_expfit_claude.py
# a_t dependence of the AXIAL tp const+exp fit, ell=1,2,3, on the ensembles present at BOTH a_t:
#   L1 g1.0  and  L2 g2.0   (Nf 2/4/6).
# SAME PHYSICAL-TIME fit window for both spacings (NM 2026-08-18): t = [2.0, 5.0] R
#   = the L1/L2 axial tp window dt[10,25] at a_t=0.2, converted per a_t (a_t=0.2 -> dt[10,25];
#     a_t=0.1 -> dt[20,50]).  m_eff and m0 are PHYSICAL (acosh/a_t), so the two spacings are
#     directly comparable; agreement => temporal continuum under control.
# const+exp fit m_eff(dt)=m0+A exp(-c1*dt_phys); config-jackknife on m0; kmin=20.  One png per ell.
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
KMIN = 20
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
PAIRS = [(1, 1.0, "1.000000"), (2, 2.0, "1.000000")]
ATS = [0.2, 0.1]
atsty = {0.2: ("o", "-", "a_t=0.2"), 0.1: ("s", "--", "a_t=0.1")}
WIN_PHYS = (2.0, 5.0)   # physical fit window t=dt*a_t (R units), shared across a_t and ell
ELLS = [1, 2, 3]


def esn(nf, L, g, at, hb):
    return ("data_Nf%d_gsq%.6fat%.6fnu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, at, L, hb))


def conn_files(nf, L, g, at, hb):
    fs = sorted(glob.glob(esn(nf, L, g, at, hb) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if int(re.search(r"corr\.(\d+)\.", f).group(1)) >= KMIN]


def load(fs, key):
    o = []
    for fn in fs:
        with h5py.File(fn, 'r') as f:
            o.append(f[key + '/real'][()] + 1j * f[key + '/imag'][()])
    return np.array(o)


def C_axial(fs, l):
    return 2.0 * (sum(load(fs, "h0/ylm_axial/s3/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)).real


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


def effmass_jk(samp, at):
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at   # PHYSICAL mass
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def model(x, m0, A, c1):
    return m0 + A * np.exp(-c1 * x)


def m0fit(me, mean, err, D, xphys):
    # D = dt indices; xphys = physical time at those dt (fit in physical time so c1 is comparable)
    if np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
        return None, None, None
    p0 = [mean[D[-1]], mean[D[0]] - mean[D[-1]], 1.0]
    bnd = ([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf])
    try:
        popt, _ = curve_fit(model, xphys, mean[D], sigma=err[D], p0=p0, bounds=bnd, maxfev=40000)
    except (RuntimeError, ValueError):
        return None, None, None
    H = me.shape[0]
    m0b = np.full(H, np.nan)
    for b in range(H):
        yb = me[b, D]
        if np.any(~np.isfinite(yb)):
            continue
        try:
            pb, _ = curve_fit(model, xphys, yb, sigma=err[D], p0=popt, bounds=bnd, maxfev=40000)
            m0b[b] = pb[0]
        except (RuntimeError, ValueError):
            continue
    good = np.isfinite(m0b)
    if good.sum() < H // 2:
        return None, None, popt
    mb = m0b[good]
    n = mb.size
    return mb.mean(), math.sqrt((n - 1.0) / n * np.sum((mb - mb.mean()) ** 2)), popt


dtp = np.arange(1, Nt // 2)
lo, hi = WIN_PHYS
allrows = []   # (ell, L, g, nf, at, m0, sig, ncfg)
for ell in ELLS:
    fig, axs = plt.subplots(1, len(PAIRS), figsize=(7 * len(PAIRS), 5))
    for ax, (L, g, hb) in zip(np.atleast_1d(axs).ravel(), PAIRS):
        for nf in NFS:
            for at in ATS:
                fs = conn_files(nf, L, g, at, hb)
                if len(fs) < 4:
                    continue
                me, mean, err = effmass_jk(C_axial(fs, ell).astype(float), at)
                D = dtp[(dtp * at >= lo - 1e-9) & (dtp * at <= hi + 1e-9)]
                xphys = D * at
                M, sig, popt = m0fit(me, mean, err, D, xphys)
                mk, ls, lab = atsty[at]
                tphys = dtp * at
                sel = tphys <= 6.0
                ax.errorbar(tphys[sel], mean[dtp][sel], yerr=err[dtp][sel], fmt=mk, ms=3.2,
                            capsize=1.5, color=nfcol[nf], ls="none", alpha=0.75,
                            label="Nf%d %s%s" % (nf, lab, "" if M is None else " m0=%.2f(%.2f)" % (M, sig)))
                if popt is not None:
                    xc = np.linspace(lo, hi, 60)
                    ax.plot(xc, model(xc, *popt), color=nfcol[nf], lw=1.5, ls=ls)
                if M is not None:
                    allrows.append((ell, L, g, nf, at, M, sig, len(fs)))
        ax.axvspan(lo, hi, color="gray", alpha=0.1)
        ax.set_title(r"axial tp $\ell=%d$  L%d g%.1f : $a_t$ dep (phys win t[%.1f,%.1f])" % (ell, L, g, lo, hi))
        ax.set_xlabel(r"physical time $t=dt\cdot a_t$ (R)")
        ax.set_ylabel(r"physical $m_\mathrm{eff}$ (1/R)")
        ax.set_xlim(0, 6)
        ax.set_ylim(0, {1: 12, 2: 16, 3: 20}[ell])
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6, ncol=1)
    fig.suptitle(r"AXIAL tp $\ell=%d$: $a_t$ dependence, const+exp on shared physical window" % ell)
    fig.tight_layout()
    fig.savefig("figs/effmass_axial_tp_at_l%d_claude.png" % ell, dpi=150)
    print("# wrote effmass_axial_tp_at_l%d_claude.png" % ell)

# ---- table: physical m0, at=0.2 vs at=0.1, per (ell,L,g,nf) ----
out = ["# Axial tp a_t dependence: physical m0 (const+exp, shared physical window t[%.1f,%.1f])"
       % (lo, hi), "",
       "| ell | L | gsq | Nf | m0(at=0.2) | m0(at=0.1) | ratio 0.1/0.2 |",
       "|-----|---|-----|----|-----------|-----------|---------------|"]
by = {}
for ell, L, g, nf, at, M, sig, n in allrows:
    by[(ell, L, g, nf, at)] = (M, sig)
seen = sorted({(ell, L, g, nf) for (ell, L, g, nf, at) in by})
print("\n# ell L g Nf : m0(0.2)  m0(0.1)  ratio")
for ell, L, g, nf in seen:
    a2 = by.get((ell, L, g, nf, 0.2))
    a1 = by.get((ell, L, g, nf, 0.1))
    s2 = "%.3f(%.3f)" % a2 if a2 else "--"
    s1 = "%.3f(%.3f)" % a1 if a1 else "--"
    r = "%.3f" % (a1[0] / a2[0]) if (a1 and a2 and a2[0]) else "--"
    print("  l%d L%d g%.1f Nf%d : %s  %s  %s" % (ell, L, g, nf, s2, s1, r))
    out.append("| %d | %d | %.1f | %d | %s | %s | %s |" % (ell, L, g, nf, s2, s1, r))
open("axial_tp_at_dependence_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote axial_tp_at_dependence_claude.md")
