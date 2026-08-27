# effmass_axial_tp_l3_fit_claude.py
# m-AVERAGED axial tp ell=3 effective mass with the const+exp fit CURVE + window shaded, L1..L4.
# This is the ell=3 correlator the l3/l1 ratio fits (ratio_axial_tp_l3l1_vs_a2_claude.py): same
# C = 2 Re<Vpp> summed over m /(2l+1), same cosh effmass, same CANONICAL windows win(3,L) =
# L1/L2 [10,25], L3/L4 [5,15].  m_eff(dt) = m0 + A exp(-c1 dt), config-jackknife, kmin=20.
# (The main effmass_axial_tp_expfit plot keeps ell=3 to L1,L2; THIS one shows L3,L4 too so the noisy
# L3 point in the l3/l1 ratio can be inspected with its fit.)  Nf 2/4/6 overlaid; one png, 4 panels.
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
LVAL = 3
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
gls = {0.5: "-", 1.0: "--", 1.5: ":", 2.0: "--", 3.0: ":", 4.5: "-.", 4.0: "--", 6.0: "-."}
FITDT_BASE = {1: (10, 25), 2: (10, 25), 3: (5, 15), 4: (5, 15)}   # canonical (no ell=3 overrides)
dtp = np.arange(1, Nt // 2)


def win(L):
    lo, hi = FITDT_BASE[L]
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
    if np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
        return None, None, None
    x = D.astype(float)
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
    return mb.mean(), math.sqrt((n - 1.0) / n * np.sum((mb - mb.mean()) ** 2)), popt


fig, axs = plt.subplots(1, 4, figsize=(24, 5))
for ax, L in zip(np.atleast_1d(axs).ravel(), [1, 2, 3, 4]):
    D = win(L)
    xc = np.linspace(D[0], D[-1], 60)
    for g in GS[L]:
        for nf in NFS:
            fs = conn_files(nf, L, g)
            if len(fs) < 4:
                continue
            me, mean, err = effmass_jk(C_axial(fs, LVAL).astype(float))
            M, sig, popt = m0fit(me, mean, err, D)
            ax.errorbar(dtp, mean[dtp], yerr=err[dtp], fmt='o', ms=2.5, capsize=1.5,
                        color=nfcol[nf], ls="none", alpha=0.6, label="Nf%d g%.1f (n=%d)" % (nf, g, len(fs)))
            if M is not None:
                ax.plot(xc, model(xc, *popt), color=nfcol[nf], lw=1.6, ls=gls[g])
                ax.axhline(M, color=nfcol[nf], ls=gls[g], lw=0.7, alpha=0.5)
    ax.axvspan(D[0], D[-1], color="gray", alpha=0.12)
    ax.set_title(r"axial tp $\ell=3$ (m-avg) L%d  fit dt[%d,%d]" % (L, D[0], D[-1]))
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 6)
    ax.set_xlabel("dt")
    ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=2)
fig.suptitle(r"AXIAL tp $\ell=3$ (m-avg): const+exp fit $m_0+A e^{-c_1 dt}$ (curve + window shaded) -- L1..L4")
fig.tight_layout()
fig.savefig("figs/effmass_axial_tp_l3_fit_claude.png", dpi=150)
print("# wrote effmass_axial_tp_l3_fit_claude.png")
