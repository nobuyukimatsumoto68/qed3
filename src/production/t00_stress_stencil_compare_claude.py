#!/usr/bin/env python3
# t00_stress_stencil_compare_claude.py
# Overlay the free-limit T_00 single-operator effmass for the O(a^2) and O(a^4) symmetric derivatives,
# plus the block-Hankel ground, against the free references sigma=2/R, T00=3/R, 2m=4/R.
# Run:  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_stress_stencil_compare_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "final/analysis_axial")
import numpy as np
import distill_contract_claude as dc
import t00_stress_claude as t0s
import effmass_axial_tp_l3_perm_hankel_claude as hk

DTMAX = 28
t0s.DTMAX = DTMAX
ROVER = 1.0 / 0.189

dual = dc.dual_areas_from_mesh()
w = dual * dc.Y00
k = dc.KS[0]

C2, _ = t0s.corr_one_config(k, w, t0s.STENCILS[2])
C4, _ = t0s.corr_one_config(k, w, t0s.STENCILS[4])


def logeff(C):
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


em2 = logeff(C2)
em4 = logeff(C4)

# block-Hankel ground on the O(a^2) correlator (positive-decaying)
emH, _ = hk.hankel_effmass_scalar(-C2, [0, 2, 4], 4, 2, 3, 0.2)
tH = np.arange(emH.shape[0])

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(9.0, 5.8))
for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R=0.378$", "tab:green"),
                      (3.0 / ROVER, r"$T_{00}=3/R=0.567$", "tab:red"),
                      (4.0 / ROVER, r"$2m=4/R=0.756$", "gray")):
    ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
    ax.text(0.3, val + 0.006, lab, fontsize=9, color=col)

dt = np.arange(DTMAX - 1)
g2 = np.isfinite(em2) & (np.arange(DTMAX - 1) >= 2)
g4 = np.isfinite(em4) & (np.arange(DTMAX - 1) >= 2)
ax.plot(dt[g2], em2[g2], color="tab:red", marker="o", ms=5, lw=1.1, label=r"single op, $O(a^2)$ deriv")
ax.plot(dt[g4], em4[g4], color="tab:orange", marker="D", ms=4.5, lw=1.1, ls="--",
        label=r"single op, $O(a^4)$ deriv")
gH = np.isfinite(emH[:, 0])
ax.plot(tH[gH], emH[gH, 0], color="tab:blue", marker="s", ms=4.5, lw=1.0,
        label=r"block-Hankel ground (off024 reb2@4 T0=3)")

ax.set_ylim(0.3, 1.0)
ax.set_xlim(2, 24)
ax.set_xlabel(r"$dt$")
ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
ax.set_title(r"Free-limit $T_{00}$ effmass: derivative order + block-Hankel  (L1 Nv=24 exact, 1 cfg)",
             fontsize=11)
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
fig.tight_layout()
os.makedirs("figs", exist_ok=True)
out = "figs/t00_stress_stencil_compare_free_claude.png"
fig.savefig(out, dpi=140)
plt.close(fig)
print("# min(O(a^2))=%.4f @ dt=%d" % (np.nanmin(em2[2:20]), 2 + int(np.nanargmin(em2[2:20]))))
print("# min(O(a^4))=%.4f @ dt=%d" % (np.nanmin(em4[2:20]), 2 + int(np.nanargmin(em4[2:20]))))
print("# Hankel ground min=%.4f @ t=%d" % (np.nanmin(emH[3:12, 0]), 3 + int(np.nanargmin(emH[3:12, 0]))))
print("# -> %s" % out)
