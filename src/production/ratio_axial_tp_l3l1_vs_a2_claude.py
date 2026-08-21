# ratio_axial_tp_l3l1_vs_a2_claude.py
# DIMENSIONLESS axial tp ratio  m0(ell=3, m-avg) / m0(ell=1, m-avg)  vs lattice spacing a^2 ~ 1/L^2.
# Self-contained CORRELATED config-jackknife: ell=1 and ell=3 are fit on the SAME delete-1 samples,
# so the ratio is formed per sample (cancels 1/a_t and captures the l3/l1 correlation).  Uses the
# CANONICAL axial tp const+exp windows (same win() as effmass_axial_tp_expfit_claude.py).
#   m_eff(dt) = m0 + A exp(-c1 dt), config-jackknife.  All L1-L4 (NM 2026-08-18: ell=3 re-included at
#   L3,L4 for this ratio even though its effmass plot is L1,L2 only).  Global ncfg>=100 cut.
# Free/CFT axial ratio Delta_(l=3)/Delta_(l=1) = 4/2 = 2  (Delta_l = l+1).
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
NCFG_MIN = 100   # global cut (NM 2026-08-18)
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
# canonical axial tp const+exp windows (dt), identical to effmass_axial_tp_expfit_claude.py.
FITDT_BASE = {1: (10, 25), 2: (10, 25), 3: (5, 15), 4: (5, 15)}
FITDT_OVR = {(1, 3): (5, 20), (1, 4): (5, 20), (2, 2): (10, 21)}
FREE = 2.0   # free/CFT axial ratio Delta_(l=3)/Delta_(l=1) = 4/2 = 2
EXCLUDE = {(4, 6.0, 6)}   # drop L4 g6.0 Nf6: err 0.31, off-trend low (NM 2026-08-18)
invL2 = {1: 1.0, 2: 0.25, 3: 1.0 / 9.0, 4: 0.0625}
dxNf = {2: -0.012, 4: 0.0, 6: 0.012}
dtp = np.arange(1, Nt // 2)


def win(ell, L):
    lo, hi = FITDT_OVR.get((ell, L), FITDT_BASE[L])
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


def meff_samples(samp):
    # per delete-1 jackknife sample: LATTICE effmass me[b] = acosh (a_t cancels in the ratio)
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)])
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def model(dt, m0, A, c1):
    return m0 + A * np.exp(-c1 * dt)


def m0_persample(me, mean, err, D):
    if np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
        return None
    x = D.astype(float)
    p0 = [mean[D[-1]], mean[D[0]] - mean[D[-1]], 0.3]
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


def ratio_jk(nf, L, g):
    fs = conn_files(nf, L, g)
    if len(fs) < NCFG_MIN:
        return None
    me1, mn1, er1 = meff_samples(C_axial(fs, 1).astype(float))
    me3, mn3, er3 = meff_samples(C_axial(fs, 3).astype(float))
    m1 = m0_persample(me1, mn1, er1, win(1, L))
    m3 = m0_persample(me3, mn3, er3, win(3, L))
    if m1 is None or m3 is None:
        return None
    good = np.isfinite(m1) & np.isfinite(m3) & (m1 > 0)
    if good.sum() < len(m1) // 2:
        return None
    Rb = m3[good] / m1[good]
    n = Rb.size
    R = Rb.mean()
    sig = math.sqrt((n - 1.0) / n * np.sum((Rb - R) ** 2))
    return R, sig, len(fs)


# ---------- collect ----------
ratios = []   # (L, g, nf, R, sig, ncfg)
for L in [1, 2, 3, 4]:
    for g in GS[L]:
        for nf in NFS:
            if (L, g, nf) in EXCLUDE:
                continue
            res = ratio_jk(nf, L, g)
            if res is None:
                continue
            R, sig, n = res
            ratios.append((L, g, nf, R, sig, n))
            print("  L%d g%.1f Nf%d : ratio l3/l1 = %.4f(%.4f)  ncfg=%d" % (L, g, nf, R, sig, n))

# ---------- plot vs 1/L^2 ----------
fig, ax = plt.subplots(figsize=(7.5, 5.5))
for nf in NFS:
    first = True
    for L in [1, 2, 3, 4]:
        for g in GS[L]:
            d = [x for x in ratios if x[0] == L and x[2] == nf and abs(x[1] - g) < 1e-9]
            if not d:
                continue
            R, sig = d[0][3], d[0][4]
            ax.errorbar(invL2[L] + dxNf[nf], R, yerr=sig, marker="o", ms=6, capsize=3,
                        color=nfcol[nf], ls="none", label=("Nf%d" % nf) if first else None)
            first = False
ax.axhline(FREE, color="gray", ls="--", lw=1.0)
ax.text(0.02, FREE, r" free $2$", color="gray", fontsize=9, va="bottom",
        transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$a^2 \sim 1/L^2$")
ax.set_ylabel(r"$m_0(\ell{=}3)/m_0(\ell{=}1)$  (axial tp, m-avg)")
ax.set_title(r"axial tp $\ell{=}3/\ell{=}1$ vs $a^2$ (const+exp $m_0$, correlated jk, ncfg$\geq$100)")
ax.set_xlim(-0.05, 1.1)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("figs/ratio_axial_tp_l3l1_vs_a2_claude.png", dpi=150)
print("# wrote ratio_axial_tp_l3l1_vs_a2_claude.png")

# ---------- md ----------
out = ["# Axial tp ratio m0(l=3)/m0(l=1) vs a^2~1/L^2 (const+exp tp windows, correlated jk, ncfg>=100)", "",
       "| L | 1/L^2 | gsq | Nf | ratio l3/l1 | ncfg |",
       "|---|-------|-----|----|-------------|------|"]
for L, g, nf, R, sig, n in sorted(ratios):
    out.append("| %d | %.4f | %.1f | %d | %.4f(%.4f) | %d |" % (L, invL2[L], g, nf, R, sig, n))
open("ratio_axial_tp_l3l1_vs_a2_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_axial_tp_l3l1_vs_a2_claude.md")
