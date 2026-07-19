# effmass_conn_claude.py
# Effective-mass analysis of the redo massless CONN tower (config-jackknife, cosh effmass).
# Group (1) AXIAL tp: l=1,2.   Group (2) SCALAR: PS & FS, l=0,1.
# One png per (group, L, gsq); Nf=2,4,6 overlaid.  kmin=20 therm cut (dir=trajectory, stride 20
# so effectively all conn configs kept; applied for safety on the config index).
import glob, math, re
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 4: "0.400000-1.000000"}
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
dt = np.arange(1, Nt // 2)
tt = dt * at


def esn(nf, L, g):
    return ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB[L]))


def files(nf, L, g):
    fs = sorted(glob.glob(esn(nf, L, g) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
    return [f for f in fs if int(re.search(r"corr\.(\d+)\.", f).group(1)) >= KMIN]


def load(fs, key):
    out = []
    for fn in fs:
        with h5py.File(fn, "r") as f:
            out.append(f[key + "/real"][()] + 1j * f[key + "/imag"][()])
    return np.array(out)


def gl(fs, grp, l, sub="Vpp"):
    return sum(load(fs, "h0/%s/l%d/m%d/%s" % (grp, l, m, sub)) for m in range(-l, l + 1)) / (2 * l + 1)


def gl_curr(fs, l, axial):
    grp = "ylm_axial" if axial else "ylm"
    return sum(load(fs, "h0/%s/s3/l%d/m%d/Vpp" % (grp, l, m)) for m in range(-l, l + 1)) / (2 * l + 1)


def C_axial(fs, l):
    return 2.0 * gl_curr(fs, l, True).real
def C_PS(fs, l):
    return 2.0 * gl(fs, "scalar", l, "Vpp").real
def C_FS(fs, l):
    return (gl(fs, "scalar", l, "Vpp") + gl(fs, "scalar_fs", l, "Vmm")).real


def meff_acosh(C):
    m = np.full(Nt, np.nan)
    for t in range(1, Nt - 1):
        den = 2.0 * C[t]
        if den == 0:
            continue
        r = (C[t - 1] + C[t + 1]) / den
        if r > 1.0:
            m[t] = np.arccosh(r)
    return m


# per-(group, L) constant-fit plateau windows in t = a_t n_t (NM 2026-07-18)
FITW = {('tp', 1): (4.0, 5.2), ('tp', 2): (2.4, 4.0),
        ('sc', 1): (3.6, 4.8), ('sc', 2): (2.4, 3.2)}


def effmass_jk(samp):
    # samp: (ncfg, Nt) real correlator -> per-sample jackknife effmass me (H,Nt), mean, err (1/a_t)
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)          # delete-1 config means
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return me, mean, err


def const_fit(me, mean, err, lo, hi):
    # inverse-variance constant fit over [lo,hi] (diagonal chi2), jackknife error via the
    # linear estimator M^b = sum_t w_t me[b,t].  Returns (M, sigma) or (None, None).
    H = me.shape[0]
    win = (tt >= lo) & (tt <= hi)
    idx = dt[win]
    good = idx[np.isfinite(mean[idx]) & (err[idx] > 0)]
    if len(good) < 2:
        return None, None
    w = 1.0 / err[good] ** 2
    w = w / w.sum()
    Mb = (me[:, good] * w).sum(1)
    M = Mb.mean()
    sig = math.sqrt((H - 1) * np.mean((Mb - M) ** 2))
    return M, sig


def panel(ax, series, title, lo, hi):
    # series: list of (label, color, linestyle/marker, C_func, l); returns # series drawn.
    # constant fit per series over [lo,hi], fit line drawn over that window.
    n_drawn = 0
    for lab, col, mk, cfun, l in series:
        for nf in NFS:
            fs = files(nf, L, g)
            if len(fs) < 4:
                continue
            samp = cfun(fs, l).astype(float)
            me, mean, err = effmass_jk(samp)
            ax.errorbar(tt, mean[dt], yerr=err[dt], fmt=mk, ms=3, capsize=2,
                        color=nfcol[nf], label="%s Nf%d" % (lab, nf))
            M, sig = const_fit(me, mean, err, lo, hi)
            if M is not None:
                ax.hlines(M, lo, hi, color=nfcol[nf], lw=1.5)
                ax.fill_between([lo, hi], M - sig, M + sig, color=nfcol[nf], alpha=0.15, linewidth=0)
                FITS.append((L, g, lab, nf, M, sig, lo, hi, len(fs)))
            n_drawn += 1
    ax.set_xlabel(r"$t = a_t n_t$")
    ax.set_ylabel(r"$m_\mathrm{eff}$ (lattice, $1/a_t$)")
    ax.set_title(title)
    ax.set_ylim(0, 6)
    ax.set_xlim(0, 8)
    ax.grid(alpha=0.3)
    if n_drawn:
        ax.legend(fontsize=6, ncol=3)
    return n_drawn


FITS = []
for L in [1, 2]:
    for g in GS[L]:
        # group (1): axial tp l=1,2
        tlo, thi = FITW[('tp', L)]
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        n1 = panel(ax, [("A l1", None, "o", C_axial, 1), ("A l2", None, "s", C_axial, 2)],
                   "axial tp effmass: L%d g%.1f (l=1,2; fit [%.1f,%.1f])" % (L, g, tlo, thi), tlo, thi)
        fn = "effmass_axial_L%d_g%.1f_claude.png" % (L, g)
        if n1:
            fig.tight_layout()
            fig.savefig(fn, dpi=150)
        plt.close(fig)
        # group (2): scalar PS/FS l=0,1
        slo, shi = FITW[('sc', L)]
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        n2 = panel(ax, [("PS l0", None, "o", C_PS, 0), ("PS l1", None, "s", C_PS, 1),
                        ("FS l0", None, "^", C_FS, 0), ("FS l1", None, "v", C_FS, 1)],
                   "scalar effmass: L%d g%.1f (PS/FS, l=0,1; fit [%.1f,%.1f])" % (L, g, slo, shi), slo, shi)
        fn2 = "effmass_scalar_L%d_g%.1f_claude.png" % (L, g)
        if n2:
            fig.tight_layout()
            fig.savefig(fn2, dpi=150)
        plt.close(fig)
        print("L%d g%.1f: axial=%d scalar=%d series%s" % (L, g, n1, n2, "" if (n1 or n2) else " (SKIP)"))

# write fitted constants to md
lines = ["# Effective-mass constant fits (redo massless conn, L1/L2)", "",
         "a_t m_eff = cosh effective mass, constant fit (inverse-variance, config-jackknife err),",
         "kmin=20. Axial tp l=1,2; scalar PS/FS l=0,1.  Nf=2,4,6.  Per-(group,L) fit windows",
         "in t=a_t n_t (the 'fit range' column). Value = a_t m (jackknife error in parens).", "",
         "| L | gsq | channel | Nf | a_t m | fit range | ncfg |",
         "|---|-----|---------|----|-------|-----------|------|"]
for L, g, lab, nf, M, sig, lo, hi, n in FITS:
    lines.append("| %d | %.1f | %s | %d | %.3f(%.3f) | [%.1f, %.1f] | %d |"
                 % (L, g, lab, nf, M, sig, lo, hi, n))
open("effmass_fits_claude.md", "w").write("\n".join(lines) + "\n")
print("# wrote effmass_fits_claude.md (%d fits)" % len(FITS))
