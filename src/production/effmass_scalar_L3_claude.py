# effmass_scalar_L3_claude.py
# Scalar PS/FS effmass at L=3 for ell=0 and ell=1, to inspect the noisy strong-coupling ensembles.
# Grid: rows = g^2 (1.5,3,4.5), cols = ell (0,1).  PS (circle) and FS (square) physical effmass vs dt,
# Nf2/4/6 overlaid.  Fit window dt[6,16] shaded (the window used in scalar_over_axial).
import glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
L = 3
HB = "0.400000-1.000000"
GS = [1.5, 3.0, 4.5]
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
WINL = {0: (6, 20), 1: (6, 13)}   # L3 windows: ell=0 [6,20], ell=1 [6,13] (NM 2026-08-18)
dtp = np.arange(1, Nt // 2)


def conn_files(nf, g):
    esn = ("data_Nf%d_gsq%.6fat0.200000nu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
           "_vmRe0.000000vmIm0.000000/" % (nf, g, L, HB))
    fs = sorted(glob.glob(esn + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
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
    return mean, err


chmk = {"PS": "o", "FS": "s"}
fig, axs = plt.subplots(len(GS), 2, figsize=(13, 4.2 * len(GS)), sharex=True)
for ir, g in enumerate(GS):
    for ic, ell in enumerate((0, 1)):
        ax = axs[ir, ic]
        for nf in NFS:
            fs = conn_files(nf, g)
            if len(fs) < 4:
                continue
            for chan in ("PS", "FS"):
                if ell == 0 and chan == "FS":
                    continue   # PS0==FS0, skip the duplicate
                mean, err = effmass_jk(C_scalar(fs, ell, chan).astype(float))
                ax.errorbar(dtp, mean[dtp], err[dtp], marker=chmk[chan], ms=3, capsize=1.5, lw=0.7,
                            color=nfcol[nf], mfc="none" if chan == "FS" else nfcol[nf],
                            label="%s Nf%d (n=%d)" % (chan, nf, len(fs)))
        w = WINL[ell]
        ax.axvspan(w[0], w[1], color="gray", alpha=0.12)
        ax.set_title(r"scalar $\ell{=}%d$  L3 $g^2$=%.1f%s" % (ell, g, "  (PS0=FS0)" if ell == 0 else ""))
        ax.set_xlim(0, 26)
        ax.set_ylim(0, 5)
        if ir == len(GS) - 1:
            ax.set_xlabel("dt")
        ax.set_ylabel(r"$m_\mathrm{eff}$ (physical)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6, ncol=2)
fig.suptitle(r"SCALAR PS/FS effmass at L3 ($\ell{=}0$ left, $\ell{=}1$ right; window dt[6,16] shaded)")
fig.tight_layout()
fig.savefig("figs/effmass_scalar_L3_claude.png", dpi=150)
print("# wrote effmass_scalar_L3_claude.png")
