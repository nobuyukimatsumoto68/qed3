# effmass_axial_l3_expfit_claude.py
# Const+exp effmass fit  m_eff(dt) = m0 + A exp(-c1 dt)  of the AXIAL-current l=3 (m-averaged, conn)
# cosh effmass.  m0 = plateau mass.  Per-L fit window (dt = timeslice sep): L1 [14,22], L2 [8,16].
# Config-jackknife (fit each delete-1 sample -> m0^b).  kmin=20.  Nf 2/4/6 overlaid; L1|L2 panels.
import glob, math, re
import numpy as np, h5py
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0]}
HB = {1: "1.000000", 2: "1.000000"}
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
FITDT = {1: np.arange(8, 21), 2: np.arange(6, 14)}   # per-L const+exp window (dt)
dtp = np.arange(1, Nt // 2)
LVAL = 3


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


# axial l=1 const-fit plateau (t = a_t n_t): established windows
AXW1 = {1: (4.0, 5.2), 2: (2.4, 4.0)}
tt = dtp * at


def axial_l1(nf, L, g):
    fs = conn_files(nf, L, g)
    if len(fs) < 4:
        return None, None
    me, mean, err = effmass_jk(C_axial(fs, 1).astype(float))
    lo, hi = AXW1[L]
    idx = dtp[(tt >= lo) & (tt <= hi)]
    good = idx[np.isfinite(mean[idx]) & (err[idx] > 0)]
    if len(good) < 2:
        return None, None
    w = 1.0 / err[good] ** 2
    w /= w.sum()
    Mb = (me[:, good] * w).sum(1)
    M = Mb.mean()
    H = me.shape[0]
    sig = math.sqrt((H - 1.0) / H * np.sum((Mb - M) ** 2))
    return M, sig


# ---------- axial l=3 const+exp fit + effmass plot ----------
gls = {0.5: "-", 1.0: "--", 1.5: ":", 2.0: "--", 3.0: ":"}
rows3 = []   # (L, g, nf, m0, sig)
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
for ax, L in zip(axs, [1, 2]):
    D = FITDT[L]
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
                ax.plot(xc, model(xc, *popt), color=nfcol[nf], lw=1.4, ls=gls[g])
                ax.axhline(M, color=nfcol[nf], ls=gls[g], lw=0.7, alpha=0.5)
                rows3.append((L, g, nf, M, sig, len(fs)))
    ax.axvspan(D[0], D[-1], color="gray", alpha=0.1)
    ax.set_title(r"axial $\ell=3$ (m-avg) L%d  fit dt[%d,%d]" % (L, D[0], D[-1]))
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 6)
    ax.set_xlabel("dt")
    ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=3)
fig.suptitle(r"AXIAL $\ell=3$: fit $m_0+A\,e^{-c_1 dt}$  (m0 = plateau)")
fig.tight_layout()
fig.savefig("effmass_axial_l3_expfit_claude.png", dpi=150)
print("# wrote effmass_axial_l3_expfit_claude.png")

lines = ["# Axial l=3 const+exp fit masses (m-avg conn, m0+A exp(-c1 dt), kmin=20)", "",
         "a_t m0 from m_eff(dt)=m0+A exp(-c1 dt), config-jackknife. windows: L1 dt[14,22], L2 dt[8,16].", "",
         "| L | gsq | Nf | a_t m0 | ncfg |", "|---|-----|----|--------|------|"]
print("\n# L gsq Nf : a_t m0 (err)  ncfg")
for L, g, nf, M, sig, n in rows3:
    print("  L%d g%.1f Nf%d : %.4f(%.4f)  %d" % (L, g, nf, M, sig, n))
    lines.append("| %d | %.1f | %d | %.4f(%.4f) | %d |" % (L, g, nf, M, sig, n))
open("axial_l3_expfit_masses_claude.md", "w").write("\n".join(lines) + "\n")
print("# wrote axial_l3_expfit_masses_claude.md (%d fits)" % len(rows3))

# ---------- ratio axial l=3 / l=1 ----------
m3 = {(L, g, nf): (M, sig) for (L, g, nf, M, sig, n) in rows3}
ratios = []   # (L, g, nf, r, er)
for L in [1, 2]:
    for g in GS[L]:
        for nf in NFS:
            if (L, g, nf) not in m3:
                continue
            m3v, m3e = m3[(L, g, nf)]
            m1v, m1e = axial_l1(nf, L, g)
            if m1v is None:
                continue
            r = m3v / m1v
            er = r * math.sqrt((m3e / m3v) ** 2 + (m1e / m1v) ** 2)
            ratios.append((L, g, nf, r, er))

FREE = 2.0   # fermion Delta_l = l+1 -> l=3/l=1 = 4/2
Lmk = {1: "o", 2: "s"}
invL2 = {1: 1.0, 2: 0.25}
dxNf = {2: -0.012, 4: 0.0, 6: 0.012}

# vs gsq
fig, ax = plt.subplots(figsize=(7, 5))
for L in [1, 2]:
    for nf in NFS:
        d = [x for x in ratios if x[0] == L and x[2] == nf]
        if not d:
            continue
        ax.errorbar([x[1] for x in d], [x[3] for x in d], [x[4] for x in d],
                    marker=Lmk[L], color=nfcol[nf], capsize=3, ls="-", lw=0.8,
                    mfc="none" if L == 2 else nfcol[nf], label="Nf%d L%d" % (nf, L))
ax.axhline(FREE, color="gray", ls="--", lw=1.0)
ax.text(0.02, FREE, r" free $2$", color="gray", fontsize=9, va="bottom", transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$g_0^2$")
ax.set_ylabel(r"$\Delta_{A,\ell=3}/\Delta_{A,\ell=1}$")
ax.set_title(r"axial $\ell=3/\ell=1$ vs $g_0^2$")
ax.grid(alpha=0.3)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig("ratio_axial_l3l1_vs_gsq_claude.png", dpi=150)
print("# wrote ratio_axial_l3l1_vs_gsq_claude.png")

# vs 1/L^2
fig, ax = plt.subplots(figsize=(7, 5))
for nf in NFS:
    for g in GS[1] + GS[2]:
        d = [x for x in ratios if x[2] == nf and abs(x[1] - g) < 1e-9]
        if not d:
            continue
        ax.errorbar([invL2[x[0]] + dxNf[nf] for x in d], [x[3] for x in d], [x[4] for x in d],
                    marker="o", color=nfcol[nf], capsize=3, ls="none",
                    label="Nf%d" % nf if abs(g - GS[1][0]) < 1e-9 else None)
ax.axhline(FREE, color="gray", ls="--", lw=1.0)
ax.text(0.02, FREE, r" free $2$", color="gray", fontsize=9, va="bottom", transform=ax.get_yaxis_transform())
ax.set_xlabel(r"$1/L^2$")
ax.set_ylabel(r"$\Delta_{A,\ell=3}/\Delta_{A,\ell=1}$")
ax.set_title(r"axial $\ell=3/\ell=1$ vs $1/L^2$")
ax.set_xlim(-0.05, 1.1)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("ratio_axial_l3l1_vs_invL2_claude.png", dpi=150)
print("# wrote ratio_axial_l3l1_vs_invL2_claude.png")
