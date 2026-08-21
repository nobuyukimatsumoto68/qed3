# effmass_axial_tp_l3_perm_claude.py
# m-DEPENDENT (per-m) axial tp ell=3 effective mass -- the individual m = -3..3 components each,
# NOT m-averaged.  Reveals any splitting among the 7 ell=3 states (on the icosahedral sphere ell=3
# reduces to T2 + G, so the m-components need not be degenerate).  Derived from the m-averaged
# effmass_axial_tp_expfit_claude.py (same C = 2 Re<Vpp>, cosh effmass, config-jackknife, kmin=20).
#   m_eff(dt) = acosh((C[t-1]+C[t+1])/(2C[t]))/a_t .  One png per L (panels = gsq), 7 m-curves + the
#   m-average overlaid (thick black).  Nf2 shown (change NF).  L1, L2 (where ell=3 is usable).
import glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
LVAL = 3
NF = 2                                 # representative flavor (per-m spread is Nf-insensitive)
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
LS = [1, 2, 3, 4]
dtp = np.arange(1, Nt // 2)
# distinct color+marker per m (color-blind: never color alone)
MS_ = list(range(-LVAL, LVAL + 1))
mcol = {-3: "tab:red", -2: "tab:orange", -1: "tab:olive", 0: "tab:green",
        1: "tab:cyan", 2: "tab:blue", 3: "tab:purple"}
mmk = {-3: "o", -2: "s", -1: "^", 0: "D", 1: "v", 2: "P", 3: "X"}


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


def C_axial_m(fs, l, m):
    # single-m axial tp correlator (2 Re<Vpp>), same normalization convention as the m-avg (sans 1/(2l+1))
    return 2.0 * load(fs, "h0/ylm_axial/s3/l%d/m%d/Vpp" % (l, m)).real


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


for L in LS:
    fig, axs = plt.subplots(1, len(GS[L]), figsize=(6.5 * len(GS[L]), 5), sharey=True)
    for ax, g in zip(np.atleast_1d(axs).ravel(), GS[L]):
        fs = conn_files(NF, L, g)
        if len(fs) < 4:
            ax.set_title(r"L%d g%.1f Nf%d: <4 cfg" % (L, g, NF))
            continue
        # per-m
        Csum = np.zeros((len(fs), Nt))
        for m in MS_:
            Cm = C_axial_m(fs, LVAL, m).astype(float)
            Csum += Cm
            mean, err = effmass_jk(Cm)
            ax.errorbar(dtp, mean[dtp], yerr=err[dtp], marker=mmk[m], ms=3, capsize=1.2, lw=0.7,
                        color=mcol[m], alpha=0.8, label="m=%+d" % m)
        # m-average (thick black) for reference
        mean, err = effmass_jk(Csum / (2 * LVAL + 1))
        ax.errorbar(dtp, mean[dtp], yerr=err[dtp], marker=".", ms=4, capsize=1.5, lw=1.4,
                    color="k", label="m-avg")
        ax.set_title(r"axial tp $\ell=3$ per-m  L%d g%.1f Nf%d (n=%d)" % (L, g, NF, len(fs)))
        ax.set_xlim(0, 30)
        ax.set_ylim(0, 6)
        ax.set_xlabel("dt")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    np.atleast_1d(axs).ravel()[0].set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    fig.suptitle(r"AXIAL tp $\ell=3$ m-DEPENDENT effmass (individual m; $\ell=3 \to T_2 + G$ on the icos sphere)  L%d" % L)
    fig.tight_layout()
    fig.savefig("figs/effmass_axial_tp_l3_perm_L%d_claude.png" % L, dpi=150)
    print("# wrote effmass_axial_tp_l3_perm_L%d_claude.png" % L)
