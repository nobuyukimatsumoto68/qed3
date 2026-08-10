# therm_traj_check_claude.py
# Flow-trajectory analysis for the redo ensembles (NM 2026-07-18 simplification):
#   fig 1: E_s(t) vs t (plain trajectory; log y -- E_s drops ~100x over the window)
#   fig 2: t^2 E_s(t) vs t with the dashed line at c; crossing = t0, g_R^2 = c/sqrt(t0)
# t = lab flow time (eps=0.01, save every step, tmax=2.0; frozen protocol).
# Data: production h5 (therm_flow_claude.h5) if present, else the v6 check .dat.
# Late-config mean (NLATE) as the equilibrated estimate; jackknife comes with the
# production analysis.
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ENS = [
    ("Nf2_gsq0.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000", "L1 g0.5"),
    ("Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000", "L1 g1.0"),
    ("Nf2_gsq1.500000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000", "L1 g1.5"),
    ("Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "L2 g1.0"),
    ("Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "L2 g2.0"),
    ("Nf2_gsq3.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000", "L2 g3.0"),
    ("Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000", "L4 g2.0"),
    ("Nf2_gsq4.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000", "L4 g4.0"),
    ("Nf2_gsq6.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000", "L4 g6.0"),
]
NLATE = 5
C = 0.003

colors = ["tab:red", "tab:orange", "tab:brown",
          "tab:blue", "tab:cyan", "tab:purple",
          "tab:green", "tab:olive", "black"]
# curves only; linestyle varies with coupling rank within each L color family (color-blind rule)
lstyles = ["-", "--", ":",
           "-", "--", ":",
           "-", "--", ":"]


def read_traj(ens):
    # production h5 preferred; v6 check .dat fallback. Returns (tlist, E[cfg][t]).
    h5path = "data_%s/therm_flow_claude.h5" % ens
    if os.path.exists(h5path):
        import h5py
        f = h5py.File(h5path, "r")
        tlist = f["tlist"][()]
        E = f["E"][()].T   # stored [t][cfg] -> [cfg][t]
        return tlist, E
    path = "data_%s/therm_series_v6_claude.dat" % ens
    tmax = None
    with open(path) as f:
        for line in f:
            if line.startswith("# tmax"):
                tmax = float(line.split("=")[1])
                break
    dat = np.loadtxt(path)
    E = dat[:, 1:]
    tlist = np.linspace(0.0, tmax, E.shape[1])
    return tlist, E


def crossing(tl, y, c):
    # first upward crossing of y(t) = c, linear interpolation
    for j in range(1, len(y)):
        if y[j - 1] < c <= y[j]:
            frac = (c - y[j - 1]) / (y[j] - y[j - 1])
            return tl[j - 1] + frac * (tl[j] - tl[j - 1])
    return None


fig, ax = plt.subplots(figsize=(7, 5))
fig2, ax2 = plt.subplots(figsize=(7, 5))
print("# late-%d-config mean; c = %g" % (NLATE, C))
print("# %-9s %8s %10s" % ("ens", "t0", "gR^2"))
for i, (ens, lab) in enumerate(ENS):
    tlist, E = read_traj(ens)
    Elate = E[-NLATE:, :].mean(axis=0)
    ax.plot(tlist, Elate, linestyle=lstyles[i], color=colors[i], label=lab, lw=1.5)
    ax2.plot(tlist, tlist**2 * Elate, linestyle=lstyles[i], color=colors[i], label=lab, lw=1.5)
    t0 = crossing(tlist, tlist**2 * Elate, C)
    if t0 is None:
        print("# %-9s NO crossing below tmax" % lab)
    else:
        print("# %-9s %8.4f %10.5f" % (lab, t0, C / np.sqrt(t0)))

ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$E_s(t)$")
ax.set_yscale("log")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("figs/therm_traj_E_claude.png", dpi=150)
print("# wrote therm_traj_E_claude.png")

ax2.set_xlabel(r"$t$")
ax2.set_ylabel(r"$t^{2}\, E_s(t)$")
ax2.axhline(C, color="gray", linestyle="--", lw=1.0, label=r"$c=%g$" % C)
ax2.set_xlim(0.0, 1.0)
ax2.set_ylim(0.0, 0.035)
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)
fig2.tight_layout()
fig2.savefig("figs/therm_traj_t2E_claude.png", dpi=150)
print("# wrote therm_traj_t2E_claude.png")
