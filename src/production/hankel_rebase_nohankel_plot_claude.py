#!/usr/bin/env python3
# hankel_rebase_nohankel_plot_claude.py
#   Plain GEVP (NO block-Hankel, NO rebase) for the 0++ two-sigma / two-meson channel,
#   drawn in the SAME scale / grid style as the shortlist figures, with the 2m_{PS} line.
#   "No Hankel" == single time-block: offsets=[0].  A full-rank rebase (nkeep = n_op) is just a
#   basis rotation, so the generalized eigenvalues (masses) equal the plain 2-op GEVP.
#   Reuses the cached store + engine from hankel_rebase_scan_claude.py (imported, not modified).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import hankel_rebase_scan_claude as hs

store, twin, ncfg, tag = hs.build_store()
M2PS = 0.644     # 2 m_{PS}, m_{PS} ~ 0.322
T0 = 3
os.makedirs("figs", exist_ok=True)

# no Hankel: single block [0]; full-keep 2-op "rebase" = plain 2x2 GEVP (masses unchanged)
em_c, em_err, ems, nbin = hs.eval_config(store, twin, 1, [0], [(5, 2)], T0)
tmax = em_c.shape[0]
nk = em_c.shape[1]
ts = np.arange(tmax)

cols = ["tab:green", "tab:red", "tab:blue"]
mkr = ["o", "s", "^"]
labs = ["state 0 (one-meson-rich ground)", "state 1 (two-meson)", "state 2"]

fig, ax = plt.subplots(figsize=(8.6, 5.6))
ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.7)
ax.text(tmax * 0.62, M2PS + 0.012, r"$2m_{PS}\approx0.644$", fontsize=9, color="gray")
for n in range(nk):
    g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.2)
    ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 3], marker=mkr[n % 3],
                ms=5, lw=1.1, capsize=2.5, label=labs[n] if n < 3 else "state %d" % n)
ax.set_ylim(0.2, 0.9)
ax.set_xlim(T0, min(tmax, 22))
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
ax.set_title("NO HANKEL  offsets=[0]  2-op  T0=3  (plain GEVP, no rebase)", fontsize=11)
ax.text(0.02, 0.02, "single time-block; state 1 (two-meson) plateaus late\ncompare vs Hankel picks (fast, flat m1 at 2mPS)",
        transform=ax.transAxes, fontsize=8, va="bottom",
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
fig.tight_layout()
out = "figs/hankel_scan_NOHANKEL_2op_T0-3_claude.png"
fig.savefig(out, dpi=130)
plt.close(fig)
print("# -> %s" % out)

print("\n#  t |   m0(err)        m1(err)")
for t in range(T0, min(tmax, 22)):
    print("#  %2d | %6.4f(%.4f)  %6.4f(%.4f)" % (t, em_c[t, 0], em_err[t, 0], em_c[t, 1], em_err[t, 1]))
