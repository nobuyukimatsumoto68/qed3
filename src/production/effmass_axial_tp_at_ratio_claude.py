# effmass_axial_tp_at_ratio_claude.py
# a_t dependence of the DIMENSIONLESS axial tp ratio  m0(ell=2)/m0(ell=1).
# The ratio cancels the 1/a_t normalization (and the lattice-vs-physical coupling ambiguity), so it is
# the clean temporal-discretization test: if R(a_t=0.1) == R(a_t=0.2), the spectrum ratio is a_t-robust.
# Ensembles present at both a_t: L1 g1.0, L2 g2.0 (Nf 2/4/6).
# const+exp fit m_eff(dt)=m0+A exp(-c1*t_phys) over the SHARED physical window t=[2.0,5.0]R.
# Correlated config-jackknife: ell=1 and ell=2 fit on the SAME delete-1 samples -> ratio per sample.
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
WIN_PHYS = (2.0, 5.0)
FREE = 1.5   # free/CFT axial ratio Delta_(l=2)/Delta_(l=1) = 3/2 (Delta_l = l+1)
dtp = np.arange(1, Nt // 2)


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


def meff_samples(samp, at):
    # per delete-1 jackknife sample: physical effmass me[b] (acosh/at), plus full-sample mean/err
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def model(x, m0, A, c1):
    return m0 + A * np.exp(-c1 * x)


def m0_persample(fs, l, at):
    # returns the per-jackknife-sample m0 array (nan where a sample fit fails)
    me, mean, err = meff_samples(C_axial(fs, l).astype(float), at)
    D = dtp[(dtp * at >= WIN_PHYS[0] - 1e-9) & (dtp * at <= WIN_PHYS[1] + 1e-9)]
    x = D * at
    if np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
        return None
    p0 = [mean[D[-1]], mean[D[0]] - mean[D[-1]], 1.0]
    bnd = ([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf])
    try:
        popt, _ = curve_fit(model, x, mean[D], sigma=err[D], p0=p0, bounds=bnd, maxfev=40000)
    except (RuntimeError, ValueError):
        return None
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
    return m0b


def ratio_jk(fs, at):
    m1 = m0_persample(fs, 1, at)
    m2 = m0_persample(fs, 2, at)
    if m1 is None or m2 is None:
        return None, None
    good = np.isfinite(m1) & np.isfinite(m2) & (m1 > 0)
    if good.sum() < len(m1) // 2:
        return None, None
    Rb = m2[good] / m1[good]
    n = Rb.size
    R = Rb.mean()
    sig = math.sqrt((n - 1.0) / n * np.sum((Rb - R) ** 2))
    return R, sig


rows = []   # (L, g, nf, at, R, sig, ncfg)
fig, axs = plt.subplots(1, len(PAIRS), figsize=(7 * len(PAIRS), 5), sharey=True)
for ax, (L, g, hb) in zip(np.atleast_1d(axs).ravel(), PAIRS):
    for nf in NFS:
        xs, ys, es = [], [], []
        for at in ATS:
            fs = conn_files(nf, L, g, at, hb)
            if len(fs) < 4:
                continue
            R, sig = ratio_jk(fs, at)
            if R is None:
                continue
            xs.append(at)
            ys.append(R)
            es.append(sig)
            rows.append((L, g, nf, at, R, sig, len(fs)))
        if xs:
            ax.errorbar(xs, ys, yerr=es, marker="o", ms=7, capsize=3, lw=1.2,
                        color=nfcol[nf], label="Nf%d" % nf)
    ax.axhline(FREE, color="k", ls=":", lw=1, alpha=0.6, label=r"free $3/2$")
    ax.set_title(r"axial tp $m_0(\ell{=}2)/m_0(\ell{=}1)$  L%d g%.1f" % (L, g))
    ax.set_xlabel(r"$a_t$")
    ax.set_xlim(0.05, 0.25)
    ax.set_ylim(1.0, 2.2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
np.atleast_1d(axs).ravel()[0].set_ylabel(r"$m_0(\ell{=}2)/m_0(\ell{=}1)$")
fig.suptitle(r"AXIAL tp ratio $\ell{=}2/\ell{=}1$: $a_t$ dependence (dimensionless; shared phys window t[2,5])")
fig.tight_layout()
fig.savefig("figs/axial_tp_ratio_l2l1_at_claude.png", dpi=150)
print("# wrote axial_tp_ratio_l2l1_at_claude.png")

out = ["# Axial tp ratio m0(l=2)/m0(l=1) vs a_t (dimensionless; shared phys window t[2.0,5.0])", "",
       "| L | gsq | Nf | a_t | ratio l2/l1 | ncfg |", "|---|-----|----|-----|-------------|------|"]
print("\n# L g Nf at : ratio l2/l1  ncfg")
for L, g, nf, at, R, sig, n in rows:
    print("  L%d g%.1f Nf%d at%.1f : %.4f(%.4f)  %d" % (L, g, nf, at, R, sig, n))
    out.append("| %d | %.1f | %d | %.1f | %.4f(%.4f) | %d |" % (L, g, nf, at, R, sig, n))
open("axial_tp_ratio_l2l1_at_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote axial_tp_ratio_l2l1_at_claude.md")
