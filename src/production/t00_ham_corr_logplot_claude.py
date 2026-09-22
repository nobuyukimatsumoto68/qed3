#!/usr/bin/env python3
# t00_ham_corr_logplot_claude.py
# Semilog plot of the free-limit T_00 (naive e.sigma, R=0) connected correlator C(dt), to check whether the
# large-dt effmass fall is NUMERICAL PRECISION LOSS (C flattens onto a noise floor) rather than a physical
# lighter (sigma) state (C would bend onto a shallower slope -0.378).  Reference slopes: 3/R=0.567, 2/R=0.378.
# Run:  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_ham_corr_logplot_claude.py

import os
for kk, vv in [("ENS", "free"), ("NVDIR", "distill_Nv24"), ("LREF", "1")]:
    os.environ.setdefault(kk, vv)
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import t00_stress_ham_claude as th
import t00_wilson_kernel_claude as wk

th.DTMAX = 30
ROVER = 1.0 / 0.189
W, _, _ = wk.build_W("../../geometry/data/", 1, r=0.0)      # naive e.sigma energy vertex
C = th.corr_one_config(dc.KS[0], W)
sgn = np.sign(C[2])
Cpos = sgn * C
dt = np.arange(th.DTMAX)

print("#  dt        C(dt)")
for d in range(1, th.DTMAX):
    if np.isfinite(Cpos[d]):
        print("  %2d   % .8e" % (d, Cpos[d]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
g = np.isfinite(Cpos) & (Cpos > 0)
fig, ax = plt.subplots(figsize=(8.8, 5.8))
ax.semilogy(dt[g], Cpos[g], "o-", color="tab:red", ms=5, lw=1.1, label=r"$C_{T_{00}}(dt)$ (naive $e\cdot\sigma$)")
# reference exponentials anchored at dt=8 (mid, still clean)
a = 8
E3 = 3.0 / ROVER
E2 = 2.0 / ROVER
ref3 = Cpos[a] * np.exp(-E3 * (dt - a))
ref2 = Cpos[a] * np.exp(-E2 * (dt - a))
ax.semilogy(dt, ref3, "--", color="k", lw=1.0, label=r"slope $3/R=0.567$")
ax.semilogy(dt, ref2, ":", color="tab:green", lw=1.2, label=r"slope $2/R=0.378$ ($\sigma$)")
ax.set_ylim(1e-9, 2.0)
ax.set_xlim(0, th.DTMAX - 1)
ax.set_xlabel(r"$dt$")
ax.set_ylabel(r"$C_{T_{00}}(dt)$")
ax.set_title(r"Free $T_{00}$ correlator (naive $e\cdot\sigma$, 1 cfg) -- log scale", fontsize=11)
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3, which="both")
fig.tight_layout()
out = "figs/t00_ham_corr_logplot_free_claude.png"
fig.savefig(out, dpi=140)
plt.close(fig)
print("\n# -> %s" % out)
