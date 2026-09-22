#!/usr/bin/env python3
# t00_stress_ham_hankel_claude.py
# Block-Hankel + rebase GEVP (frozen core) on the naive e.sigma (R=0) T_00 correlator, to get the free plateau
# EARLIER than the point-to-point effmass.  Correlator capped at DTMAX below the ~1e-8 precision floor.
# Params (sigma/distillation thread): offsets [0,2,4], reb2@4 (NKEEP=2 @ REBT=4), T0=3.  Single free config.
# Run:  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_stress_ham_hankel_claude.py

import os
for kk, vv in [("ENS", "free"), ("NVDIR", "distill_Nv24"), ("LREF", "1")]:
    os.environ.setdefault(kk, vv)
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "final/analysis_axial")
import numpy as np
import distill_contract_claude as dc
import t00_stress_ham_claude as th
import t00_wilson_kernel_claude as wk
import effmass_axial_tp_l3_perm_hankel_claude as hk

DTMAX = int(os.environ.get("DTMAX", "24"))          # stay below the ~1e-8 precision floor (C(24)~7e-8)
# Best plateau for the naive e.sigma T00 (NM 2026-09-16): off0-3 reb1@3 T0=2 (flat ~0.562 over t=9-15).
OFFS = [int(x) for x in os.environ.get("OFFSETS", "0,3").split(",")]
REBT = int(os.environ.get("REBT", "3"))
NKEEP = int(os.environ.get("NKEEP", "1"))
T0 = int(os.environ.get("T0", "2"))
ROVER = 1.0 / 0.189

th.DTMAX = DTMAX
W, _, _ = wk.build_W("../../geometry/data/", 1, r=0.0)
C = th.corr_one_config(dc.KS[0], W)
Cpos = np.sign(C[2]) * C

em, Vop = hk.hankel_effmass_scalar(Cpos, OFFS, REBT, NKEEP, T0, 0.2)

# point-to-point for overlay
with np.errstate(all="ignore"):
    emp = np.log(Cpos[:-1] / Cpos[1:])

tmax1 = em.shape[0]
print("# block-Hankel off%s reb%d@%d T0=%d on naive e.sigma T00 (free). refs sigma=0.378 T00=0.567"
      % ("-".join(map(str, OFFS)), NKEEP, REBT, T0))
print("#  t   ground     1st-exc    point2point")
for t in range(tmax1):
    if np.isfinite(em[t, 0]):
        pp = emp[t] if t < emp.shape[0] else np.nan
        e1 = em[t, 1] if NKEEP > 1 else np.nan
        print("  %2d  %8.4f  %8.4f   %8.4f" % (t, em[t, 0], e1, pp))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
tt = np.arange(tmax1)
fig, ax = plt.subplots(figsize=(8.8, 5.6))
for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R$", "tab:green"),
                      (3.0 / ROVER, r"$T_{00}=3/R=0.567$", "tab:red"),
                      (4.0 / ROVER, r"$2m=4/R$", "gray")):
    ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
    ax.text(0.2, val + 0.006, lab, fontsize=9, color=col)
tp = np.arange(emp.shape[0])
gp = np.isfinite(emp)
ax.plot(tp[gp], emp[gp], "-", color="lightgray", lw=1.0, marker=".", ms=4, label="point-to-point")
g0 = np.isfinite(em[:, 0])
ax.plot(tt[g0], em[g0, 0], "o-", color="tab:red", ms=5, lw=1.2, label="block-Hankel ground")
if NKEEP > 1:
    g1 = np.isfinite(em[:, 1])
    ax.plot(tt[g1], em[g1, 1], "s:", color="tab:blue", ms=4.5, lw=1.0, markerfacecolor="none",
            label="block-Hankel 1st-exc")
ax.set_ylim(0.3, 1.0)
ax.set_xlim(0, tmax1)
ax.set_xlabel(r"$dt$")
ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
ax.set_title(r"Free $T_{00}$ (naive $e\cdot\sigma$): block-Hankel plateau  off%s reb%d@%d T0=%d  1 cfg"
             % ("-".join(map(str, OFFS)), NKEEP, REBT, T0), fontsize=10.5)
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
fig.tight_layout()
out = "figs/t00_stress_ham_hankel_free_claude.png"
fig.savefig(out, dpi=140)
plt.close(fig)
print("\n# -> %s" % out)
