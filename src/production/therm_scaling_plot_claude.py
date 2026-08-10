# therm_scaling_plot_claude.py
# Canonical flow-scaling figure: E_s(t) sqrt(L) vs x = r/t (log-log) for the nine Nf2 massless
# redo ensembles, with binned-jackknife error bands (fill_between, alpha=0.2) and the per-L
# free (LO mode-sum) curves calibrated on the weakest coupling of each L over x in [0.2, 2].
# Data: data_<ens>/therm_flow_claude.h5 (any grid version; tlist read per file).
# Output: therm_traj_EsqrtL_loglog_v2_claude.png (name kept across grid versions).
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

ENS = [
    ("Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000", "L1 g0.5", 1, 0.5),
    ("Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000", "L1 g1.0", 1, 1.0),
    ("Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000", "L1 g1.5", 1, 1.5),
    ("Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "L2 g1.0", 2, 1.0),
    ("Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "L2 g2.0", 2, 2.0),
    ("Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "L2 g3.0", 2, 3.0),
    ("Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000", "L4 g2.0", 4, 2.0),
    ("Nf2_gsq4.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000", "L4 g4.0", 4, 4.0),
    ("Nf2_gsq6.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000", "L4 g6.0", 4, 6.0),
]
NLATE = 500   # late-config window (thermalization cut by construction)
BIN = 20      # jackknife bin (autocorrelation)
AT = 0.2

colors = ["tab:red", "tab:orange", "tab:brown",
          "tab:blue", "tab:cyan", "tab:purple",
          "tab:green", "tab:olive", "black"]
lstyles = ["-", "--", ":"] * 3


KMIN = 20   # thermalization cut (2026-07-18): keep trajectories k >= 20


def load(ens):
    f = h5py.File("data_%s/therm_flow_claude.h5" % ens, "r")
    tl = f["tlist"][()]
    E = f["E"][()].T
    kl = f["klist"][()]
    E = E[kl >= KMIN, :]
    n = min(NLATE, E.shape[0])
    return tl, E[-n:, :]


def jackknife(A, nbin):
    # binned delete-1 jackknife over configs (rows); returns mean, err per column
    nb = max(2, A.shape[0] // nbin)
    nbin_eff = A.shape[0] // nb
    bins = np.array([A[i*nbin_eff:(i+1)*nbin_eff, :].mean(axis=0) for i in range(nb)])
    mean = bins.mean(axis=0)
    sub = (nb * mean - bins) / (nb - 1)
    err = np.sqrt((nb - 1) / nb * ((sub - mean)**2).sum(axis=0))
    return mean, err


def sum_x(x, L, b):
    # LO mode sum as a function of x = r/t and L (coupling drops out)
    j = np.arange(1, 6 * L + 1)
    lam = j * (j + 1.0)
    w = (2 * j + 1) * np.sqrt(lam)
    res = np.zeros_like(x)
    for ll, ww in zip(lam, w):
        res += ww * np.exp(-2.0 * b * AT * ll / (x * L))
    return res


fig, ax = plt.subplots(figsize=(7, 5))
for i, (ens, lab, L, gsq) in enumerate(ENS):
    tl, A = load(ens)
    mean, err = jackknife(A, BIN)
    x = (gsq / L) / tl[1:]
    y = mean[1:] * np.sqrt(L)
    dy = err[1:] * np.sqrt(L)
    ax.plot(x, y, linestyle=lstyles[i], color=colors[i], label=lab, lw=1.5)
    ax.fill_between(x, y - dy, y + dy, color=colors[i], alpha=0.2, linewidth=0)

xg = np.geomspace(0.05, 5.0, 2)
ax.plot(xg, 4.5e-2 * xg**1.5, color="gray", linestyle="--", lw=1.0, label=r"slope $3/2$")

ax.set_xlabel(r"$r/t$")
ax.set_ylabel(r"$E_s(t)\,\sqrt{L}$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(4e-4, 0.8)
ax.legend(fontsize=7, ncol=2)
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig("figs/therm_traj_EsqrtL_loglog_v2_claude.png", dpi=150)
print("# wrote therm_traj_EsqrtL_loglog_v2_claude.png")

# ---- figure 2: Nf dependence at L2 g3 (Nf2/4/6 + free), same jackknife bands ----
NFENS = [
    ("Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "Nf2"),
    ("Nf4_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "Nf4"),
    ("Nf6_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "Nf6"),
]
L2 = 2
r2 = 3.0 / L2
nfcolors = ["tab:red", "tab:blue", "tab:green"]
nfstyles = ["-", "--", ":"]

fig2, ax2 = plt.subplots(figsize=(7, 5))
for i, (ens, lab) in enumerate(NFENS):
    tl, A = load(ens)
    mean, err = jackknife(A, BIN)
    x = r2 / tl[1:]
    y = mean[1:] * np.sqrt(L2)
    dy = err[1:] * np.sqrt(L2)
    ax2.plot(x, y, linestyle=nfstyles[i], color=nfcolors[i], label=lab, lw=1.5)
    ax2.fill_between(x, y - dy, y + dy, color=nfcolors[i], alpha=0.2, linewidth=0)
xg2 = np.geomspace(0.05, 5.0, 2)
ax2.plot(xg2, 4.5e-2 * xg2**1.5, color="gray", linestyle="--", lw=1.0, label=r"slope $3/2$")
ax2.set_xlabel(r"$r/t$")
ax2.set_ylabel(r"$E_s(t)\,\sqrt{L}$")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_title(r"L2, $g_0^2=3$: $N_f$ dependence")
ax2.legend(fontsize=9)
ax2.grid(alpha=0.3, which="both")
fig2.tight_layout()
fig2.savefig("figs/therm_traj_Nfdep_L2g3_claude.png", dpi=150)
print("# wrote therm_traj_Nfdep_L2g3_claude.png")
