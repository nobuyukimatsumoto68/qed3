# ratio_vs_at2_claude.py
# a_t dependence: masses of (ell=2 tp axial) and (F glueball, l=1) IN UNITS OF (ell=1 tp axial),
# plotted against a_t^2.  Paired half-a_t ensembles: L1 g1.0, L2 g2.0 (Nf2/4/6) at a_t=0.2 and 0.1.
#   R_A = m0(axial tp, l=2) / m0(axial tp, l=1)   (both fermionic; correlated config-jackknife)
#   R_F = m0(F, l=1)        / m0(axial tp, l=1)    (F from GEVP fit_perm; naive error propagation)
# Axial: cosh effmass, const+exp m0 over the SHARED PHYSICAL window t=[2,5]R (so the two a_t are
#   compared consistently).  F: fit_perm variance-avg over t=[0.2,1.0] (physical; same at both a_t).
#   a_t=0.2 gevp tag "Nf_g_L"; a_t=0.1 tag "Nf_g_L_at0.1" (from run_glue_gevp_at01_claude.sh).
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
PAIRS = [(1, 1.0), (2, 2.0)]        # (L, gsq) paired ensembles
ATS = [0.2, 0.1]
# axial shared physical window, PER ell: ell=1 [2,5]; ell=2 [1,2.4] (NM 2026-08-18: the ell=2 late
# tail is noise-floor at a_t=0.1, so fit its earlier plateau).
# axial physical window, PER a_t (NM 2026-08-18): at=0.1 -> [1,2.4] (its late tail is noise-floor),
# at=0.2 -> [1.4,4.0].  Same window for ell=1 and ell=2.
WIN_PHYS_AT = {0.1: (1.0, 2.4), 0.2: (1.4, 4.0)}
WIN_OVR = {(0.1, 2): (1.2, 2.6)}    # (at,ell) override: at0.1 ell=2 -> [1.2,2.6] (NM 2026-08-18)
FWIN = (0.2, 1.0)                   # F l=1 fit_perm window (physical, central)
dtp = np.arange(1, Nt // 2)


def esn(nf, L, g, at):
    hb = "1.000000"
    suf = "_vmRe0.000000vmIm0.000000" if at == 0.2 else ""   # at0.1 dirs: hb1.000000 (no vm suffix)
    # try both variants (some at0.2 dirs carry the vm suffix, some not)
    return ("data_Nf%d_gsq%.6fat%.6fnu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            % (nf, g, at, L, hb))


def conn_files(nf, L, g, at):
    base = esn(nf, L, g, at)
    fs = []
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


def meff_samples(samp, at):
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at   # physical effmass
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def model(x, m0, A, c1):
    return m0 + A * np.exp(-c1 * x)


def m0_persample(fs, l, at):
    me, mean, err = meff_samples(C_axial(fs, l).astype(float), at)
    wlo, whi = WIN_OVR.get((at, l), WIN_PHYS_AT[at])
    D = dtp[(dtp * at >= wlo - 1e-9) & (dtp * at <= whi + 1e-9)]
    x = D * at
    if len(D) < 3 or np.any(~np.isfinite(mean[D])) or np.any(err[D] <= 0):
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


def _mean_err(m0b):
    g = m0b[np.isfinite(m0b) & (m0b > 0)]
    if g.size < len(m0b) // 2:
        return None
    m = g.mean()
    return m, math.sqrt((g.size - 1.0) / g.size * np.sum((g - m) ** 2))


def axial_masses(nf, L, g, at):
    # returns (m1, e1, m2, e2) PHYSICAL axial l=1,l=2 masses (None if a fit fails)
    fs = conn_files(nf, L, g, at)
    if len(fs) < 4:
        return None
    r1 = m0_persample(fs, 1, at)
    r2 = m0_persample(fs, 2, at)
    a1 = _mean_err(r1) if r1 is not None else None
    a2 = _mean_err(r2) if r2 is not None else None
    return a1, a2


def glue_mass(pref, win, nf, L, g, at):
    tag = "Nf%d_g%.6f_L%d%s" % (nf, g, L, "" if at == 0.2 else "_at0.1")
    jk = "%s_%s_jk_claude.dat" % (pref, tag)
    try:
        open(jk).close()
    except OSError:
        return None, None
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(win[0]), str(win[1]), "0,1,2"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "variance-avg" in line:
            seg = line.split("M =")[-1] if "==>" in line else line.split("M=")[-1]
            m = seg.strip().split()[0]
            try:
                return float(m.split("(")[0]), float(m.split("(")[1].rstrip(")"))
            except (ValueError, IndexError):
                return None, None
    return None, None


def F_mass(nf, L, g, at):
    return glue_mass("gevp_F", FWIN, nf, L, g, at)          # F l=1, [0.2,1.0]


def F_l2_mass(nf, L, g, at):
    return glue_mass("gevp_Fl2", (0.2, 0.6), nf, L, g, at)  # F l=2 (H), [0.2,0.6]


# ---------- collect ----------
rows = []   # (L, g, nf, at, RA, eRA, RF, eRF)
for L, g in PAIRS:
    for nf in NFS:
        for at in ATS:
            am = axial_masses(nf, L, g, at)
            a1, a2 = (None, None) if am is None else am
            mF = F_mass(nf, L, g, at)
            mF = mF if (mF is not None and mF[0] is not None) else None
            mF2 = F_l2_mass(nf, L, g, at)
            mF2 = mF2 if (mF2 is not None and mF2[0] is not None) else None
            # BARE masses in lattice a_t units (m*a_t = acosh): m_phys * a_t
            b1 = (a1[0] * at, a1[1] * at) if a1 else None
            b2 = (a2[0] * at, a2[1] * at) if a2 else None
            bF = (mF[0] * at, mF[1] * at) if mF else None
            bF2 = (mF2[0] * at, mF2[1] * at) if mF2 else None
            rows.append((L, g, nf, at, b1, b2, bF, bF2))
            print("  L%d g%.1f Nf%d at%.1f : m*at  axl1=%s axl2=%s Fl1=%s Fl2=%s" %
                  (L, g, nf, at,
                   ("%.4f(%.4f)" % b1) if b1 else "--",
                   ("%.4f(%.4f)" % b2) if b2 else "--",
                   ("%.4f(%.4f)" % bF) if bF else "--",
                   ("%.4f(%.4f)" % bF2) if bF2 else "--"))

# ---------- plot: 2 panels (paired ensembles), BARE masses m*a_t vs a_t ----------
SER = [(4, "o", "-", True, "axial l1"), (5, "^", "--", True, "axial l2"),
       (6, "s", ":", False, "F l1"), (7, "D", "-.", False, "F l2")]
fig, axs = plt.subplots(1, len(PAIRS), figsize=(7 * len(PAIRS), 5.5), sharey=True)
for ax, (L, g) in zip(np.atleast_1d(axs).ravel(), PAIRS):
    for nf in NFS:
        d = sorted([r for r in rows if r[0] == L and r[2] == nf], key=lambda r: r[3])
        for idx, mk, ls, fill, lab in SER:
            dd = [r for r in d if r[idx] is not None]
            if not dd:
                continue
            ax.errorbar([r[3] for r in dd], [r[idx][0] for r in dd], [r[idx][1] for r in dd],
                        marker=mk, ms=7, capsize=3, color=nfcol[nf], ls=ls, lw=1.0,
                        mfc=nfcol[nf] if fill else "none", label="Nf%d %s" % (nf, lab))
    ax.set_title(r"L%d  $g^2$=%.1f" % (L, g))
    ax.set_xlabel(r"$a_t$")
    ax.set_xlim(0.0, 0.25)
    ax.set_ylim(0.1, 0.55)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=3)
np.atleast_1d(axs).ravel()[0].set_ylabel(r"$m\,a_t$  (lattice $\mathrm{acosh}$)")
fig.suptitle(r"$a_t$ dependence: BARE masses $m\,a_t$ of axial tp $\ell{=}1,2$ and F $\ell{=}1$ vs $a_t$  "
             r"(axial $a_t$-flat, F $\propto a_t$)")
fig.tight_layout()
fig.savefig("figs/ratio_vs_at2_claude.png", dpi=150)
print("# wrote ratio_vs_at2_claude.png")

# ---------- SEPARATE log-log plot: m*a_t vs a_t (slope 1 = m*at prop a_t = glue; slope 0 = flat) ----------
fig, axs = plt.subplots(1, len(PAIRS), figsize=(7 * len(PAIRS), 5.5), sharey=True)
for ax, (L, g) in zip(np.atleast_1d(axs).ravel(), PAIRS):
    for nf in NFS:
        d = sorted([r for r in rows if r[0] == L and r[2] == nf], key=lambda r: r[3])
        for idx, mk, ls, fill, lab in SER:
            dd = [r for r in d if r[idx] is not None]
            if not dd:
                continue
            ax.errorbar([r[3] for r in dd], [r[idx][0] for r in dd], [r[idx][1] for r in dd],
                        marker=mk, ms=7, capsize=3, color=nfcol[nf], ls=ls, lw=1.0,
                        mfc=nfcol[nf] if fill else "none", label="Nf%d %s" % (nf, lab))
    # slope guides: prop a_t (glue) and flat (fermion), anchored at a_t=0.2
    ax.plot([0.1, 0.2], [0.2, 0.4], color="gray", ls=":", lw=1.0, alpha=0.6)   # slope 1 (m*at ~ a_t)
    ax.plot([0.1, 0.2], [0.35, 0.35], color="gray", ls=":", lw=1.0, alpha=0.6)  # slope 0 (flat)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(r"L%d  $g^2$=%.1f  (log-log)" % (L, g))
    ax.set_xlabel(r"$a_t$")
    ax.set_xlim(0.08, 0.25)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, ncol=3)
np.atleast_1d(axs).ravel()[0].set_ylabel(r"$m\,a_t$  (lattice $\mathrm{acosh}$)")
fig.suptitle(r"$a_t$ dependence (LOG-LOG): $m\,a_t$ of axial tp $\ell{=}1,2$ and F $\ell{=}1,2$ vs $a_t$  "
             r"(guides: slope 1 $=$ glue $\propto a_t$; slope 0 $=$ fermion flat)")
fig.tight_layout()
fig.savefig("figs/ratio_vs_at2_loglog_claude.png", dpi=150)
print("# wrote ratio_vs_at2_loglog_claude.png")

# ---------- md ----------
out = ["# a_t dependence: BARE masses m*a_t (lattice acosh) of axial l1,l2 & F l1,l2 vs a_t (paired ens)", "",
       "m*a_t = m_phys * a_t.  axial const+exp (at0.1 win [1,2.4], at0.2 [1.4,4.0]; at0.1 ell2 [1.2,2.6]);",
       "F l1 fit_perm [0.2,1.0], F l2 [0.2,0.6].", "",
       "| L | gsq | Nf | a_t | axl1*at | axl2*at | Fl1*at | Fl2*at |",
       "|---|-----|----|-----|---------|---------|--------|--------|"]
for L, g, nf, at, b1, b2, bF, bF2 in sorted(rows):
    s1 = "%.4f(%.4f)" % b1 if b1 else "--"
    s2 = "%.4f(%.4f)" % b2 if b2 else "--"
    sF = "%.4f(%.4f)" % bF if bF else "--"
    sF2 = "%.4f(%.4f)" % bF2 if bF2 else "--"
    out.append("| %d | %.1f | %d | %.1f | %s | %s | %s | %s |" % (L, g, nf, at, s1, s2, sF, sF2))
open("ratio_vs_at2_claude.md", "w").write("\n".join(out) + "\n")
print("# wrote ratio_vs_at2_claude.md")
