#!/usr/bin/env python3
# hankel_rebase_dt036_fit_claude.py
#   Correlated constant fit of the [0,3,6] 2-op T0=3 rebase 2@t=4 effmass plateaus.
#   state 0: t in [4,6] ; state 1 (two-meson): t in [12,14].
#   Correlated fit uses the jackknife covariance over the window:
#     mhat = (1^T Cinv y)/(1^T Cinv 1),  var = 1/(1^T Cinv 1),  chi2 = (y-mhat)^T Cinv (y-mhat).
#   Same figure style; fit bands overlaid.

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
OFF = [0, 3, 6]
WIN = {0: (4, 6), 1: (12, 14)}

em_c, em_err, ems, nbin = hs.eval_config(store, twin, 1, OFF, [(REBT, 2)], T0)
tmax = em_c.shape[0]


def corr_const_fit(lev, tlo, thi):
    ts = np.arange(tlo, thi + 1)
    y = em_c[ts, lev]
    Y = ems[:, ts, lev]                                  # (nbin, npts)
    Ym = Y.mean(0)
    # jackknife covariance of the mean
    Cov = (nbin - 1.0) / nbin * np.einsum("ip,iq->pq", Y - Ym, Y - Ym)
    Cinv = np.linalg.inv(Cov)
    one = np.ones(len(ts))
    denom = one @ Cinv @ one
    mhat = (one @ Cinv @ y) / denom
    var = 1.0 / denom
    r = y - mhat
    chi2 = r @ Cinv @ r
    dof = len(ts) - 1
    # uncorrelated (diagonal) cross-check
    w = 1.0 / em_err[ts, lev] ** 2
    m_unc = (w * y).sum() / w.sum()
    e_unc = 1.0 / np.sqrt(w.sum())
    return mhat, np.sqrt(var), chi2, dof, m_unc, e_unc


print("# [0,3,6] 2-op T0=3 rebase 2@t=4   nbin=%d (binsize 10)" % nbin)
print("# level | window |   correlated fit        chi2/dof |  uncorr (diag) fit")
fits = {}
for lev, (tlo, thi) in WIN.items():
    mhat, err, chi2, dof, m_unc, e_unc = corr_const_fit(lev, tlo, thi)
    fits[lev] = (mhat, err, tlo, thi)
    print("#   %d   | [%2d,%2d] |  %7.4f +- %6.4f   %6.2f/%d |  %7.4f +- %6.4f"
          % (lev, tlo, thi, mhat, err, chi2, dof, m_unc, e_unc))
print("# 2m_PS ~ %.4f ; m0/m1 = %.3f" % (M2PS, fits[0][0] / fits[1][0]))

# ---- plot with fit bands ----
ts = np.arange(tmax)
cols = ["tab:green", "tab:red", "tab:blue"]
mkr = ["o", "s", "^"]
labs = ["state 0 (one-meson-rich ground)", "state 1 (two-meson)"]
fig, ax = plt.subplots(figsize=(8.6, 5.6))
ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.7)
ax.text(tmax * 0.62, M2PS + 0.012, r"$2m_{PS}\approx0.644$", fontsize=9, color="gray")
for n in range(2):
    g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.2)
    ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n], marker=mkr[n],
                ms=5, lw=1.1, capsize=2.5, label=labs[n])
    mhat, err, tlo, thi = fits[n]
    ax.fill_between([tlo - 0.3, thi + 0.3], mhat - err, mhat + err, color=cols[n], alpha=0.25)
    ax.plot([tlo - 0.3, thi + 0.3], [mhat, mhat], color=cols[n], lw=1.6)
    ax.text(thi + 0.4, mhat + (0.02 if n == 0 else -0.03),
            r"$%.4f(%d)$" % (mhat, round(err * 1e4)), color=cols[n], fontsize=9)
ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.4)
ax.set_ylim(0.2, 0.9)
ax.set_xlim(T0, min(tmax, 22))
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
ax.set_title("Dt=3,6  [0,3,6]  2-op  T0=3  rebase 2@t=4  (constant fits)", fontsize=11)
ax.text(0.02, 0.02,
        "state0 [4,6]: %.4f(%d)\nstate1 [12,14]: %.4f(%d)" %
        (fits[0][0], round(fits[0][1] * 1e4), fits[1][0], round(fits[1][1] * 1e4)),
        transform=ax.transAxes, fontsize=9, va="bottom",
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
fig.tight_layout()
out = "figs/hankel_scan_Dt036_2op_T0-3_reb4_fit_claude.png"
fig.savefig(out, dpi=130)
plt.close(fig)
print("# -> %s" % out)
