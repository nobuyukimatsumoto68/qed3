# effmass_axial_tp_l3_perm_traj_claude.py
# Per-m axial tp ell=3 effmass along the a=0.5 CONSTANT-PHYSICS trajectory (g^2 = 0.5*L), L=1..4 side by
# side -- so the m-dependent ell=3 structure can be tracked toward the continuum (a^2 ~ 1/L^2 -> 0).
#   trajectory: L1 g0.5, L2 g1.0, L3 g1.5, L4 g2.0 (Nf2).  Same per-m effmass as
#   effmass_axial_tp_l3_perm_claude.py: individual m=-3..3 (color+marker) + m-average (thick black).
import glob, math, re
import numpy as np, h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nt = 128
at = 0.2
KMIN = 20
LVAL = 3
NF = 2
GTRAJ = {1: 0.5, 2: 1.0, 3: 1.5, 4: 2.0}     # a=0.5 trajectory: g^2 = 0.5*L
HB = {1: "1.000000", 2: "1.000000", 3: "0.400000-1.000000", 4: "0.400000-1.000000"}
LS = [1, 2, 3, 4]
invL2 = {1: 1.0, 2: 0.25, 3: 1.0 / 9.0, 4: 0.0625}
dtp = np.arange(1, Nt // 2)
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


fig, axs = plt.subplots(1, len(LS), figsize=(5.5 * len(LS), 5), sharey=True)
for ax, L in zip(np.atleast_1d(axs).ravel(), LS):
    g = GTRAJ[L]
    fs = conn_files(NF, L, g)
    if len(fs) < 4:
        ax.set_title(r"L%d g%.1f: <4 cfg" % (L, g))
        continue
    Csum = np.zeros((len(fs), Nt))
    for m in MS_:
        Cm = C_axial_m(fs, LVAL, m).astype(float)
        Csum += Cm
        mean, err = effmass_jk(Cm)
        ax.errorbar(dtp, mean[dtp], yerr=err[dtp], marker=mmk[m], ms=3, capsize=1.2, lw=0.7,
                    color=mcol[m], alpha=0.8, label="m=%+d" % m)
    mean, err = effmass_jk(Csum / (2 * LVAL + 1))
    ax.errorbar(dtp, mean[dtp], yerr=err[dtp], marker=".", ms=4, capsize=1.5, lw=1.4,
                color="k", label="m-avg")
    ax.set_title(r"L%d  $g^2$=%.1f  $a^2$=%.3f (n=%d)" % (L, g, invL2[L], len(fs)))
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 6)
    ax.set_xlabel("dt")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
np.atleast_1d(axs).ravel()[0].set_ylabel(r"$a_t\, m_\mathrm{eff}$")
fig.suptitle(r"AXIAL tp $\ell=3$ m-dependent effmass along the $a{=}0.5$ trajectory ($g^2{=}0.5L$), Nf%d -- L1..L4 (coarse $\to$ fine)" % NF)
fig.tight_layout()
fig.savefig("figs/effmass_axial_tp_l3_perm_traj_a05_claude.png", dpi=150)
print("# wrote effmass_axial_tp_l3_perm_traj_a05_claude.png")
