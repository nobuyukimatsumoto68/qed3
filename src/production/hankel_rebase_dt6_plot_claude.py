#!/usr/bin/env python3
# hankel_rebase_dt6_plot_claude.py
#   Dt=1,2,3,4,5,6 block-Hankel (offsets [0..6]), 2-op base, T0=3, rebase 2@t=4.
#   Same scale / grid style as the shortlist figures, with the interacting 2m_{PS} line.

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
M2PS = 0.644
T0 = 3
REBT = 4
OFF = [0, 1, 2, 3, 4, 5, 6]
os.makedirs("figs", exist_ok=True)

em_c, em_err, ems, nbin = hs.eval_config(store, twin, 1, OFF, [(REBT, 2)], T0)
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
ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.5)
ax.set_ylim(0.2, 0.9)
ax.set_xlim(T0, min(tmax, 22))
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
ax.set_title("Dt=1..6  offsets=[0..6]  2-op  T0=3  rebase 2@t=4", fontsize=11)
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
fig.tight_layout()
out = "figs/hankel_scan_Dt123456_2op_T0-3_reb4_claude.png"
fig.savefig(out, dpi=130)
plt.close(fig)
print("# -> %s  (tmax=%d)" % (out, tmax))

print("\n#  t |   m0(err)        m1(err)")
for t in range(T0, min(tmax, 22)):
    print("#  %2d | %6.4f(%.4f)  %6.4f(%.4f)" % (t, em_c[t, 0], em_err[t, 0], em_c[t, 1], em_err[t, 1]))
