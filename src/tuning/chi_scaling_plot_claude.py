#!/usr/bin/env python3
# chi_scaling_plot_claude.py
# _claude: volume-scaling SSB test (chi_volume_scaling_ssb_claude.md).  Two figures:
#   (1) chi = (1/V) sum_low 1/|z|^2  vs V     -- flat = symmetric, rising ~V = SSB.
#   (2) lambda_min vs V (log-log)             -- slope -1/2 (1/L) = symmetric, -1 (1/V) = SSB.
# Reads the sig1 (wall M=1) Arnoldi dats; overlap eigenvalue z = 1 + (mu-1)/|mu-1| (M=1).
# One representative reasonable-gsq point per L (config mean +/- std).  Run: python3 chi_scaling_plot_claude.py

import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nlow = 20

def load(path):
    re = []
    im = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        re.append(float(p[1]))
        im.append(float(p[2]))
    return np.array(re), np.array(im)

def overlap_absz(re, im):
    mu = re + 1j*im
    a = mu - 1.0
    z = 1.0 + a/np.abs(a)
    return np.abs(z)

# representative reasonable-gsq ensemble per L: (V, glob)
ENS = {
    1: (1536,  "eig_dw_arnoldi_L1_sig1_gsq1.0_cfg*_nt128_claude.dat"),
    2: (5376,  "eig_dw_arnoldi_L2_sig1_gsq1.0_cfg*_nt128_claude.dat"),
    4: (20736, "eig_dw_arnoldi_L4_sig1_gsq2.0_cfg*_nt128_claude.dat"),
}

Ls = []
Vs = []
lam_m = []
lam_e = []
chi_m = []
chi_e = []
for L in (1, 2, 4):
    V, pat = ENS[L]
    files = sorted(glob.glob(pat))
    lam = []
    chi = []
    for f in files:
        az = np.sort(overlap_absz(*load(f)))
        lam.append(az[0])
        chi.append(np.sum(1.0/az[:Nlow]**2)/V)
    lam = np.array(lam)
    chi = np.array(chi)
    Ls.append(L)
    Vs.append(V)
    lam_m.append(lam.mean())
    lam_e.append(lam.std(ddof=1)/np.sqrt(len(lam)) if len(lam) > 1 else 0.0)
    chi_m.append(chi.mean())
    chi_e.append(chi.std(ddof=1)/np.sqrt(len(chi)) if len(chi) > 1 else 0.0)

Vs = np.array(Vs, float)
lam_m = np.array(lam_m)
lam_e = np.array(lam_e)
chi_m = np.array(chi_m)
chi_e = np.array(chi_e)

# ---- Figure 1: chi vs V (the SSB fork) ----
fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.errorbar(Vs, chi_m, yerr=chi_e, marker="o", ms=8, color="C3", ls="none", capsize=4,
            label=r"$\chi=\frac{1}{V}\sum_{\rm low}1/\lambda^2$ (data)")
# SSB reference: chi ~ V, anchored at the L1 point
ssb = chi_m[0]*(Vs/Vs[0])
ax.plot(Vs, ssb, ls="--", color="C0", lw=1.2, label=r"SSB: $\chi\propto V$")
ax.axhline(chi_m[0], ls=":", color="0.4", lw=1.2, label="symmetric: flat")
ax.set_xscale("log")
ax.set_xlabel(r"$V=N_s N_t$")
ax.set_ylabel(r"$\chi$")
ax.set_ylim(0.0, max(chi_m.max()*1.4, ssb.max()*0.1))
ax.set_title(r"Low-mode chiral susceptibility vs volume (massless, reasonable gsq)")
ax.legend(loc="upper left", fontsize=9)
ax.grid(True, alpha=0.25)
fig.tight_layout()
fig.savefig("chi_vs_V_claude.png", dpi=130)
print("wrote chi_vs_V_claude.png")

# ---- Figure 2: lambda_min vs V (log-log, 1/L vs 1/V slopes) ----
fig, ax = plt.subplots(figsize=(6.4, 5.0))
ax.errorbar(Vs, lam_m, yerr=lam_e, marker="s", ms=8, color="C0", ls="none", capsize=4,
            label=r"$\lambda_{\min}$ (data)")
sym = lam_m[0]*(Vs/Vs[0])**(-0.5)
ssb = lam_m[0]*(Vs/Vs[0])**(-1.0)
ax.plot(Vs, sym, ls="--", color="C2", lw=1.2, label=r"symmetric $\propto V^{-1/2}$ ($1/L$)")
ax.plot(Vs, ssb, ls=":", color="C3", lw=1.4, label=r"SSB $\propto V^{-1}$ ($1/V$)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$V=N_s N_t$")
ax.set_ylabel(r"$\lambda_{\min}$ (overlap)")
ax.set_title(r"Overlap spectral gap vs volume")
ax.legend(loc="lower left", fontsize=9)
ax.grid(True, which="both", alpha=0.25)
fig.tight_layout()
fig.savefig("lammin_vs_V_claude.png", dpi=130)
print("wrote lammin_vs_V_claude.png")

print()
print("chi(V):", dict(zip([int(v) for v in Vs], np.round(chi_m, 4))))
print("lam(V):", dict(zip([int(v) for v in Vs], np.round(lam_m, 4))))
