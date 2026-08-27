# mat_vs_scales_claude.py
# The a_t-artifact finding: plot m*a_t (= lattice acosh, the dimensionless per-timeslice decay) for
# F(l=1) and Axial-tp(l=1) against three x-axes: a_t, g^2*R, g^2*a.  R=1, a=R/L (so g^2R=g^2, g^2a=g^2/L).
# Paired half-a_t ensembles L1 g1.0 & L2 g2.0 (Nf2/4/6), a_t in {0.2,0.1}.
# KEY: F's m*a_t ~ proportional to a_t (glue scales properly -> line through ~origin); Axial's m*a_t is
# ~a_t-INDEPENDENT (fermion artifact -> flat).  m*a_t = m_phys * a_t.
#   Axial l=1 m_phys: cosh effmass const+exp over shared PHYSICAL window t=[2,5]R.
#   F l=1 m_phys:     GEVP fit_perm variance-avg over t=[0.2,1.0] (gevp_F, at0.1 tag _at0.1).
import subprocess, glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
KMIN = 20
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
PAIRS = [(1, 1.0), (2, 2.0)]      # (L, gsq)
ATS = [0.2, 0.1]
# axial physical window per a_t (NM 2026-08-18): at=0.1 -> [1,2.4], at=0.2 -> [1.4,4.0].
WIN_PHYS_AT = {0.1: (1.0, 2.4), 0.2: (1.4, 4.0)}
FWIN = (0.2, 1.0)
dtp = np.arange(1, Nt // 2)


def conn_files(nf, L, g, at):
    base = "data_Nf%d_gsq%.6fat%.6fnu01.000000mRe0.000000mIm0.000000nt128L%d_hb1.000000" % (nf, g, at, L)
    for suf in ("", "_vmRe0.000000vmIm0.000000"):
        fs = sorted(glob.glob(base + suf + "/corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
        if fs:
            break
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


def model(x, m0, A, c1):
    return m0 + A * np.exp(-c1 * x)


def axial_phys(nf, L, g, at):
    # returns physical m0(axial l=1) with jk error, const+exp over physical window [2,5]
    fs = conn_files(nf, L, g, at)
    if len(fs) < 4:
        return None
    samp = C_axial(fs, 1).astype(float)
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    wlo, whi = WIN_PHYS_AT[at]
    D = dtp[(dtp * at >= wlo - 1e-9) & (dtp * at <= whi + 1e-9)]
    if len(D) < 3 or np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
        return None
    x = D * at
    p0 = [mean[D[-1]], mean[D[0]] - mean[D[-1]], 1.0]
    bnd = ([0.0, -np.inf, 0.0], [np.inf, np.inf, np.inf])
    try:
        popt, _ = curve_fit(model, x, mean[D], sigma=err[D], p0=p0, bounds=bnd, maxfev=40000)
    except (RuntimeError, ValueError):
        return None
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
        return None
    mb = m0b[good]
    return mb.mean(), math.sqrt((len(mb) - 1.0) / len(mb) * np.sum((mb - mb.mean()) ** 2))


def F_phys(nf, L, g, at):
    tag = "Nf%d_g%.6f_L%d%s" % (nf, g, L, "" if at == 0.2 else "_at0.1")
    jk = "gevp_F_%s_jk_claude.dat" % tag
    try:
        open(jk).close()
    except OSError:
        return None
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(FWIN[0]), str(FWIN[1]), "0,1,2"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "variance-avg" in line:
            seg = line.split("M =")[-1]
            m = seg.strip().split()[0]
            try:
                return float(m.split("(")[0]), float(m.split("(")[1].rstrip(")"))
            except (ValueError, IndexError):
                return None
    return None


# ---------- collect: y = m*a_t (lattice acosh) ----------
rows = []   # (chan, L, g, nf, at, y, ey, g2R, g2a)
for L, g in PAIRS:
    for nf in NFS:
        for at in ATS:
            g2R = g            # R=1
            g2a = g / L        # a = R/L
            fa = axial_phys(nf, L, g, at)
            if fa is not None:
                rows.append(("Axial", L, g, nf, at, fa[0] * at, fa[1] * at, g2R, g2a))
            ff = F_phys(nf, L, g, at)
            if ff is not None:
                rows.append(("F", L, g, nf, at, ff[0] * at, ff[1] * at, g2R, g2a))
for r in rows:
    print("  %-5s L%d g%.1f Nf%d at%.1f : m*at=%.4f(%.4f)" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6]))

# ---------- 3 panels: y=m*a_t vs {a_t, g^2R, g^2a} ----------
XSPEC = [("$a_t$", 4, (0.0, 0.25)), (r"$g^2 R$", 7, (0.0, 2.5)), (r"$g^2 a$", 8, (0.5, 1.5))]
fig, axs = plt.subplots(1, 3, figsize=(18, 5.5), sharey=True)
for ax, (xlab, xi, xlim) in zip(axs, XSPEC):
    for chan, mk, fill in [("Axial", "s", False), ("F", "o", True)]:
        for nf in NFS:
            d = [r for r in rows if r[0] == chan and r[3] == nf]
            if not d:
                continue
            ax.errorbar([r[xi] for r in d], [r[5] for r in d], [r[6] for r in d], marker=mk, ms=7,
                        mfc=nfcol[nf] if fill else "none", capsize=3, ls="none", color=nfcol[nf],
                        label="%s Nf%d" % (chan, nf))
    ax.set_xlabel(xlab)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, None)
    ax.grid(alpha=0.3)
    ax.set_title(r"$m\,a_t$ vs %s" % xlab)
    if xi == 4:
        ax.legend(fontsize=7, ncol=2)
axs[0].set_ylabel(r"$m\,a_t$  (lattice $\mathrm{acosh}$)")
fig.suptitle(r"$m\,a_t$ for F ($\ell{=}1$, filled) and Axial tp ($\ell{=}1$, open) vs $a_t$, $g^2R$, $g^2a$  "
             r"(F $\propto a_t$; Axial $a_t$-flat = fermion artifact)")
fig.tight_layout()
fig.savefig("figs/mat_vs_scales_claude.png", dpi=150)
print("# wrote mat_vs_scales_claude.png")

out = ["# m*a_t (lattice acosh) for F(l=1) and Axial(l=1) vs a_t / g^2R / g^2a.  R=1, a=R/L.", "",
       "| chan | L | gsq | Nf | a_t | g2R | g2a | m*a_t |", "|------|---|-----|----|-----|-----|-----|-------|"]
for chan, L, g, nf, at, y, ey, g2R, g2a in rows:
    out.append("| %s | %d | %.1f | %d | %.1f | %.2f | %.2f | %.4f(%.4f) |" % (chan, L, g, nf, at, g2R, g2a, y, ey))
open("mat_vs_scales_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote mat_vs_scales_claude.md")
