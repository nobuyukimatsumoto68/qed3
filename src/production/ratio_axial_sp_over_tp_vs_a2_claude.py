# ratio_axial_sp_over_tp_vs_a2_claude.py
# Axial SP/TP mass ratios vs lattice spacing a^2 ~ 1/L^2:
#   R1 = m0(sp, ell=1) / m0(tp, ell=1)     R2 = m0(sp, ell=2) / m0(tp, ell=1)
# sp = SO(3) spatial scalar = s1 + s2 (m-summed); tp = s3 (temporal).  Both = 2 Re<Vpp> (sp is
# uniformly negative but cosh-foldable -- verified), cosh effmass acosh/a_t.  CORRELATED config-
# jackknife: sp(l1), sp(l2), tp(l1) fit on the SAME delete-1 samples -> ratios per sample (cancels
# 1/a_t and captures the sp-tp correlation).  Canonical const+exp windows, all L1-L4, ncfg>=100.
# NB: no free-value reference lines drawn (the sp/tp continuum ratio on S2xR is not assumed here).
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
NCFG_MIN = 100
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
FITDT_BASE = {1: (10, 25), 2: (10, 25), 3: (5, 15), 4: (5, 15)}
FITDT_OVR = {(1, 3): (5, 20), (1, 4): (5, 20)}     # ell=1 (both sp and tp) at L3,L4 -> [5,20]
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


def C_tp(fs, l):
    return 2.0 * (sum(load(fs, "h0/ylm_axial/s3/l%d/m%d/Vpp" % (l, m)) for m in range(-l, l + 1)) / (2 * l + 1)).real


def C_sp(fs, l):
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


def meff_samples(samp):
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)])   # LATTICE effmass (a_t cancels in ratio)
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


SP_ELLS = [0, 1]   # sp numerator ells (NM 2026-08-18: ell=0 and ell=1; ell=2,3 not needed)


def ratios_jk(nf, L, g):
    # R_l = m0(sp, ell=l) / m0(tp, ell=1) for l in SP_ELLS, correlated over the same jk samples.
    fs = conn_files(nf, L, g)
    if len(fs) < NCFG_MIN:
        return None
    me_t1, mn_t1, er_t1 = meff_samples(C_tp(fs, 1).astype(float))
    t1 = m0_persample(me_t1, mn_t1, er_t1, win(1, L))
    if t1 is None:
        return None
    out = {}
    for l in SP_ELLS:
        me_s, mn_s, er_s = meff_samples(C_sp(fs, l).astype(float))
        s = m0_persample(me_s, mn_s, er_s, win(l, L))
        if s is None:
            out[l] = None
            continue
        good = np.isfinite(t1) & np.isfinite(s) & (t1 > 0)
        if good.sum() < len(t1) // 2:
            out[l] = None
            continue
        Rb = s[good] / t1[good]
        n = Rb.size
        out[l] = (Rb.mean(), math.sqrt((n - 1.0) / n * np.sum((Rb - Rb.mean()) ** 2)))
    return out, len(fs)


# ---------- collect ----------
rows = []   # (L, g, nf, {l: (R,e)}, ncfg)
for L in [1, 2]:                        # L=3,4 dropped for now (NM 2026-08-18)
    for g in GS[L]:
        for nf in NFS:
            res = ratios_jk(nf, L, g)
            if res is None:
                continue
            out, n = res
            rows.append((L, g, nf, out, n))
            cells = "  ".join("sp%d/tp1=%s" % (l, ("%.4f(%.4f)" % out[l] if out[l] else "--")) for l in SP_ELLS)
            print("  L%d g%.1f Nf%d : %s  ncfg=%d" % (L, g, nf, cells, n))

# ---------- plot: one panel per sp ell, R_l = m0(sp,l)/m0(tp,l=1) vs 1/L^2 ----------
fig, axs = plt.subplots(1, len(SP_ELLS), figsize=(7 * len(SP_ELLS), 5.5), sharex=True)
for ax, l in zip(np.atleast_1d(axs).ravel(), SP_ELLS):
    for nf in NFS:
        first = True
        for L in [1, 2, 3, 4]:
            for g in GS[L]:
                d = [r for r in rows if r[0] == L and r[2] == nf and abs(r[1] - g) < 1e-9]
                if not d or d[0][3].get(l) is None:
                    continue
                rr = d[0][3][l]
                ax.errorbar(invL2[L] + dxNf[nf], rr[0], yerr=rr[1], marker="o", ms=6, capsize=3,
                            color=nfcol[nf], ls="none", label=("Nf%d" % nf) if first else None)
                first = False
    ax.set_xlabel(r"$a^2 \sim 1/L^2$")
    ax.set_ylabel(r"$m_0(\mathrm{sp},\ell{=}%d)/m_0(\mathrm{tp},\ell{=}1)$" % l)
    ax.set_title(r"sp $\ell{=}%d$ / tp $\ell{=}1$" % l)
    ax.set_xlim(-0.05, 1.1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
fig.suptitle(r"AXIAL sp/tp mass ratios vs $a^2$ (const+exp $m_0$, correlated jk, ncfg$\geq$100)")
fig.tight_layout()
fig.savefig("figs/ratio_axial_sp_over_tp_vs_a2_claude.png", dpi=150)
print("# wrote ratio_axial_sp_over_tp_vs_a2_claude.png")

# ---------- md ----------
hdr = "| L | 1/L^2 | gsq | Nf | " + " | ".join("sp%d/tp1" % l for l in SP_ELLS) + " | ncfg |"
sep = "|---|-------|-----|----|" + "|".join(["---------"] * len(SP_ELLS)) + "|------|"
out = ["# Axial sp/tp mass ratios vs a^2~1/L^2 (const+exp windows, correlated jk, ncfg>=100)", "",
       "R_l = m0(sp,l)/m0(tp,l=1) for l in %s.  sp=s1+s2, tp=s3." % SP_ELLS, "", hdr, sep]
for L, g, nf, o, n in sorted(rows, key=lambda r: (r[0], r[1], r[2])):
    cs = " | ".join(("%.4f(%.4f)" % o[l] if o[l] else "--") for l in SP_ELLS)
    out.append("| %d | %.4f | %.1f | %d | %s | %d |" % (L, invL2[L], g, nf, cs, n))
open("ratio_axial_sp_over_tp_vs_a2_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_axial_sp_over_tp_vs_a2_claude.md")
