#!/usr/bin/env python3
# plot_l2_modecontact_claude.py -- comparison plot of the interacting-L2 sigma^2 leak vs the MODE_CONTACT fix,
#   built from the logged effmass tables (no re-run of the heavy L2 4-vertex contraction).
#   Reads leak_L2_modecontact0_claude.log (position-space contact = bug) and leak_L2_modecontact1_claude.log (fix).
#   Series: C_11 (single meson = m_PS), C_22 position-contact (LEAKS to m_PS), C_22 mode-contact (recovers two-meson).

import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MPS = 0.3527
M2PS = 0.7054


def parse(path):
    ts = []
    c11 = []
    e11 = []
    c22 = []
    e22 = []
    with open(path) as f:
        for line in f:
            m = re.match(r"^#\s+(\d+)\s+\|\s+([-\d.]+)\(([\d.]+)\)\s+([-\d.]+)\(([\d.]+)\)", line)
            if m:
                ts.append(int(m.group(1)))
                c11.append(float(m.group(2)))
                e11.append(float(m.group(3)))
                c22.append(float(m.group(4)))
                e22.append(float(m.group(5)))
    return np.array(ts), np.array(c11), np.array(e11), np.array(c22), np.array(e22)


t0, c11_0, e11_0, c22_0, e22_0 = parse("leak_L2_modecontact0_claude.log")
t1, c11_1, e11_1, c22_1, e22_1 = parse("leak_L2_modecontact1_claude.log")

fig, ax = plt.subplots(figsize=(9.0, 6.0))
ax.axhline(M2PS, color="gray", ls="--", lw=1.1, alpha=0.8)
ax.text(11.2, M2PS + 0.012, r"$2m_{PS}=%.4f$" % M2PS, fontsize=10, color="gray")
ax.axhline(MPS, color="firebrick", ls=":", lw=1.2, alpha=0.8)
ax.text(11.2, MPS + 0.012, r"$m_{PS}=%.4f$" % MPS, fontsize=10, color="firebrick")

# single meson C_11 (same in both) -- red filled circle
g = (e11_1 < 0.1)
ax.errorbar(t1[g], c11_1[g], yerr=e11_1[g], color="firebrick", marker="o", ms=6, lw=1.2, capsize=3,
            label=r"$C_{11}$ ($\sigma_{00}$, single meson $= m_{PS}$)")
# C_22 position-space contact (the bug) -- gray filled triangle, LEAKS to m_PS
g = (e22_0 < 0.1)
ax.errorbar(t0[g], c22_0[g], yerr=e22_0[g], color="dimgray", marker="v", ms=6, lw=1.2, capsize=3,
            label=r"$C_{22}$ ($\sigma^2_{00}$), position contact -- LEAKS to $m_{PS}$")
# C_22 mode-space contact (the fix) -- blue filled square, recovers two-meson
g = (e22_1 < 0.15)
ax.errorbar(t1[g], c22_1[g], yerr=e22_1[g], color="tab:blue", marker="s", ms=6, lw=1.3, capsize=3,
            label=r"$C_{22}$ ($\sigma^2_{00}$), MODE contact (fix) -- recovers two-meson")

ax.set_ylim(0.25, 1.05)
ax.set_xlim(2, 15)
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
ax.set_title(r"Interacting $L2$ (Nf2 $g^2{=}1.0$), truncated $N_v{=}24$: the contact fix restores $\sigma^2\to$ two-meson"
             "\n" r"80 cfg, reb4 $T_0{=}2$", fontsize=10)
ax.legend(fontsize=9, loc="lower left")
ax.grid(alpha=0.3)
fig.tight_layout()
out = "figs/sigma2_L2_modecontact_fix_compare_claude.png"
fig.savefig(out, dpi=140)
print("# -> %s" % out)
