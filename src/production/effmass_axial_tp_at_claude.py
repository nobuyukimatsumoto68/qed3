# effmass_axial_tp_at_claude.py
# a_t (temporal-spacing) DEPENDENCE check for the AXIAL tp l=1 channel.
# Compares the two ensembles that exist at BOTH a_t = 0.1 and a_t = 0.2 for the same physical coupling:
#   L1 g1.0  and  L2 g2.0   (Nf 2/4/6).
# The cosh effmass acosh(...) is the LATTICE mass a_t*m; dividing by a_t gives the PHYSICAL mass m
# (1/R units), so the two spacings are directly comparable.  Plotted vs PHYSICAL time t = dt*a_t so the
# a_t=0.1 (finer, 2x more timeslices per unit t) and a_t=0.2 curves share one horizontal axis.
# Agreement => a_t -> 0 under control; a gap => temporal-discretization artifact.
# m-averaged, conn, config-jackknife, kmin=20.
import glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
KMIN = 20
NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
LVAL = 1
# paired ensembles (L, gsq, hb) present at both a_t
PAIRS = [(1, 1.0, "1.000000"), (2, 2.0, "1.000000")]
ATS = [0.2, 0.1]
atmark = {0.2: ("o", "-", "a_t=0.2"), 0.1: ("s", "--", "a_t=0.1")}


def esn(nf, L, g, at, hb):
    return ("data_Nf%d_gsq%.6fat%.6fnu01.000000mRe0.000000mIm0.000000nt128L%d_hb%s"
            "_vmRe0.000000vmIm0.000000/" % (nf, g, at, L, hb))


def conn_files(nf, L, g, at, hb):
    fs = sorted(glob.glob(esn(nf, L, g, at, hb) + "corr_ylm_conn_t00_nhits1_s1/corr.*.h0.h5"))
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
    me = np.array([meff_acosh(jk[i]) for i in range(H)]) / at   # /at -> PHYSICAL mass (1/R)
    mean = np.nanmean(me, 0)
    err = np.sqrt(np.maximum((H - 1) * np.nanmean((me - mean) ** 2, 0), 0.0))
    return mean, err


dtp = np.arange(1, Nt // 2)
fig, axs = plt.subplots(1, len(PAIRS), figsize=(7 * len(PAIRS), 5))
for ax, (L, g, hb) in zip(np.atleast_1d(axs).ravel(), PAIRS):
    for nf in NFS:
        for at in ATS:
            fs = conn_files(nf, L, g, at, hb)
            if len(fs) < 4:
                continue
            mean, err = effmass_jk(C_axial(fs, LVAL).astype(float), at)
            tphys = dtp * at
            mk, ls, lab = atmark[at]
            sel = tphys <= 4.0
            ax.errorbar(tphys[sel], mean[dtp][sel], yerr=err[dtp][sel], fmt=mk, ms=3.5,
                        capsize=1.5, color=nfcol[nf], ls=ls, lw=0.7, alpha=0.8,
                        label="Nf%d %s (n=%d)" % (nf, lab, len(fs)))
    ax.set_title(r"axial tp $\ell=1$  L%d g%.1f : $a_t$ dependence" % (L, g))
    ax.set_xlabel(r"physical time $t = dt\cdot a_t$  (units of R)")
    ax.set_ylabel(r"physical $m_\mathrm{eff}$  (1/R)")
    ax.set_xlim(0, 4)
    ax.set_ylim(0, 12)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=2)
fig.suptitle(r"AXIAL tp $\ell=1$: $a_t=0.1$ (squares/dashed) vs $a_t=0.2$ (circles/solid), physical units")
fig.tight_layout()
fig.savefig("figs/effmass_axial_tp_at_claude.png", dpi=150)
print("# wrote effmass_axial_tp_at_claude.png")
