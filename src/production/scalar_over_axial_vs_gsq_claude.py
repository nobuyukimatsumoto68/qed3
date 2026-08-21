# scalar_over_axial_vs_gsq_claude.py
# Scalar (PS, FS) at ell=0,1 in units of axial tp ell=1, vs g^2.  L=1,2.
#   PS = 2*Vpp ;  FS = Vpp + Vmm^FS.  (NB at ell=0 the two coincide: PS0 == FS0.)
#   scalar m0: cosh effmass const+exp fit, config-jackknife (physical acosh/a_t).
#   axial l=1: const+exp m0 from axial_tp_masses_summary_claude.md.  Naive error propagation.
# Two panels: PS (ell0,ell1)/axial and FS (ell0,ell1)/axial.  No free-value line drawn.
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
LLIST = [1, 2, 3, 4]
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
Lmk = {1: "o", 2: "s", 3: "D", 4: "^"}
invL2 = {1: 1.0, 2: 0.25, 3: 1.0 / 9.0, 4: 0.0625}
SUMMARY = "axial_tp_masses_summary_claude.md"
NCFG_MIN = 100   # global ncfg cut (NM 2026-08-18): drop points with <100 configs
# scalar const+exp windows (dt): ell=1 L1/L2 [8,24] L3/L4 [6,16]; ell=0 same (borrowed; validate).
WIN = {(0, 1): (8, 24), (0, 2): (8, 24), (0, 3): (6, 20), (0, 4): (6, 16),
       (1, 1): (8, 24), (1, 2): (8, 24), (1, 3): (6, 13), (1, 4): (6, 16)}
ellmk = {0: "o", 1: "^"}
dtp = np.arange(1, Nt // 2)


def axial_summary():
    d = {}
    for line in open(SUMMARY):
        s = line.strip()
        if not s.startswith("|"):
            continue
        c = [x.strip() for x in s.strip("|").split("|")]
        if len(c) != 11 or not c[0].isdigit() or int(c[0]) != 1:
            continue
        d[(int(c[2]), float(c[4]), int(c[5]))] = (float(c[8]), float(c[9]))
    return d


axsum = axial_summary()


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


def m0fit(fs, l, chan, D):
    me, mean, err = effmass_jk(C_scalar(fs, l, chan).astype(float))
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
    good = np.isfinite(m0b)
    if good.sum() < H // 2:
        return None
    mb = m0b[good]
    return mb.mean(), math.sqrt((len(mb) - 1.0) / len(mb) * np.sum((mb - mb.mean()) ** 2))


# ---------- collect ----------
rows = []   # (chan, ell, L, g, nf, r, er)
for chan in ("PS", "FS"):
    for ell in (0, 1):
        for L in LLIST:
            for nf in NFS:
                for g in GS[L]:
                    if (L, g, nf) not in axsum:
                        continue
                    fs = conn_files(nf, L, g)
                    if len(fs) < NCFG_MIN:
                        continue
                    lo, hi = WIN[(ell, L)]
                    r0 = m0fit(fs, ell, chan, np.arange(lo, hi + 1))
                    if r0 is None:
                        continue
                    a, ea = axsum[(L, g, nf)]
                    r = r0[0] / a
                    er = r * math.sqrt((r0[1] / r0[0]) ** 2 + (ea / a) ** 2)
                    rows.append((chan, ell, L, g, nf, r, er))
                    print("  %s l%d L%d g%.1f Nf%d : /axial = %.4f(%.4f)" % (chan, ell, L, g, nf, r, er))

# ---------- vs g^2: SEPARATE png per ell (each: PS + FS panels) ----------
for ell in (0, 1):
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5), sharex=True)
    for ax, chan in zip(axs, ("PS", "FS")):
        for L in LLIST:
            for nf in NFS:
                d = sorted([x for x in rows if x[0] == chan and x[1] == ell and x[2] == L and x[4] == nf],
                           key=lambda x: x[3])
                if not d:
                    continue
                ax.errorbar([x[3] for x in d], [x[5] for x in d], [x[6] for x in d],
                            marker=Lmk[L], ms=6, capsize=3, ls="-", lw=0.8, color=nfcol[nf],
                            mfc="none" if L in (2, 4) else nfcol[nf], label="Nf%d L%d" % (nf, L))
        if ell == 0:
            ax.axhline(1.0, color="gray", ls="--", lw=1.0, alpha=0.7)   # ratio=1 ref (ell=0 only)
        ax.set_title(r"scalar %s $\ell{=}%d$ / axial $\ell{=}1$" % (chan, ell))
        ax.set_ylabel(r"$\Delta_{%s,\ell=%d}/\Delta_{A,\ell=1}$" % (chan, ell))
        ax.set_xlabel(r"$g_0^2$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6, ncol=2)
    fig.suptitle(r"SCALAR (PS, FS) $\ell{=}%d$ / axial tp $\ell{=}1$ vs $g_0^2$  (L: o s D ^ = 1 2 3 4; open=L2,L4)" % ell)
    fig.tight_layout()
    fig.savefig("figs/ratio_scalar_over_axial_l%d_vs_gsq_claude.png" % ell, dpi=150)
    print("# wrote ratio_scalar_over_axial_l%d_vs_gsq_claude.png" % ell)

# ---------- vs a^2 ~ 1/L^2: SEPARATE png per ell (all g overlaid) ----------
dxNf = {2: -0.012, 4: 0.0, 6: 0.012}
for ell in (0, 1):
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5), sharex=True)
    for ax, chan in zip(axs, ("PS", "FS")):
        for nf in NFS:
            first = True
            for L in LLIST:
                for g in GS[L]:
                    d = [x for x in rows if x[0] == chan and x[1] == ell and x[2] == L
                         and x[4] == nf and abs(x[3] - g) < 1e-9]
                    if not d:
                        continue
                    ax.errorbar(invL2[L] + dxNf[nf], d[0][5], d[0][6], marker="o", ms=6, capsize=3,
                                ls="none", color=nfcol[nf], label=("Nf%d" % nf) if first else None)
                    first = False
        if ell == 0:
            ax.axhline(1.0, color="gray", ls="--", lw=1.0, alpha=0.7)   # ratio=1 ref (ell=0 only)
        ax.set_title(r"scalar %s $\ell{=}%d$ / axial $\ell{=}1$" % (chan, ell))
        ax.set_ylabel(r"$\Delta_{%s,\ell=%d}/\Delta_{A,\ell=1}$" % (chan, ell))
        ax.set_xlim(-0.05, 1.1)
        ax.set_xlabel(r"$a^2 \sim 1/L^2$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle(r"SCALAR (PS, FS) $\ell{=}%d$ / axial tp $\ell{=}1$ vs $a^2$ (1/L^2)  (all $g^2$ overlaid)" % ell)
    fig.tight_layout()
    fig.savefig("figs/ratio_scalar_over_axial_l%d_vs_a2_claude.png" % ell, dpi=150)
    print("# wrote ratio_scalar_over_axial_l%d_vs_a2_claude.png" % ell)

out = ["# Scalar (PS,FS) l=0,1 / axial l=1 vs g^2, L=1,2 (axial from summary).  PS0==FS0.", "",
       "| chan | ell | L | gsq | Nf | /axial |", "|------|-----|---|-----|----|--------|"]
for chan, ell, L, g, nf, r, er in sorted(rows):
    out.append("| %s | %d | %d | %.1f | %d | %.4f(%.4f) |" % (chan, ell, L, g, nf, r, er))
open("ratio_scalar_over_axial_vs_gsq_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_scalar_over_axial_vs_gsq_claude.md")
