# effmass_l2nf6_at_claude.py
# Axial tp effmass for L2 g2.0 Nf6, ell=1 and ell=2, a_t=0.2 vs 0.1, PHYSICAL effmass = acosh/a_t vs
# physical t = a_t*n_t.  Shows the flat plateau + noise-floor tail that destabilizes the const+exp
# const+exp fit at a_t=0.1 (198 cfg).  Window t=[2,5] shaded.
import glob, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
KMIN = 20
nf, L, g = 4, 2, 2.0
ATS = [(0.2, "tab:orange", "o"), (0.1, "tab:purple", "s")]
WIN_AT = {0.2: (1.4, 4.0), 0.1: (1.0, 2.4)}   # at-dependent fit window (shaded per a_t)
dtp = np.arange(1, Nt // 2)


def conn_files(at):
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


def effmass_jk(samp, at):
    H = samp.shape[0]
    jk = (samp.sum(0) - samp) / (H - 1)
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return mean, err


fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
for ax, l in zip(axs, [1, 2]):
    for at, col, mk in ATS:
        fs = conn_files(at)
        if len(fs) < 4:
            continue
        mean, err = effmass_jk(C_axial(fs, l).astype(float), at)
        t = dtp * at
        sel = t <= 6.0
        ax.errorbar(t[sel], mean[dtp][sel], err[dtp][sel], marker=mk, ms=4, capsize=2, lw=0.7,
                    color=col, label="at%.1f (n=%d)" % (at, len(fs)))
        w = WIN_AT[at]
        ax.axvspan(w[0], w[1], color=col, alpha=0.10)   # this a_t's fit window
    ax.set_title(r"axial tp $\ell=%d$  L2 $g^2$=2.0 Nf%d" % (l, nf))
    ax.set_xlabel(r"physical $t = a_t n_t$")
    ax.set_ylabel(r"$m_\mathrm{eff}$ (physical, acosh/$a_t$)")
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 6)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
fig.suptitle(r"L2 $g^2$2.0 Nf%d axial effmass: $a_t$0.2 vs 0.1 (at-dep window shaded: at0.2 [1.4,4], at0.1 [1,2.4])" % nf)
fig.tight_layout()
fig.savefig("figs/effmass_l2nf4_at_claude.png", dpi=150)
print("# wrote effmass_l2nf4_at_claude.png")
