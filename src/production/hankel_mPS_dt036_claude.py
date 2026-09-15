#!/usr/bin/env python3
# hankel_mPS_dt036_claude.py
#   Single-sigma pseudoscalar m_PS extracted with the SAME block-Hankel + rebase strategy used for
#   the two-meson channel, so the threshold 2 m_PS is pinned on equal footing with m_1.
#     base   : 1 operator, the single-sigma correlator M(t) = -Tr[Phi(t) tau(t,s) Phi(s) tau(s,t)]
#     Hankel : offsets [0,3,6]  (GPOF, Aubin-Orginos 1010.0202) -> 3x3, resolves a small meson tower
#     rebase : 2 states @ t=4, metric t0=3 ; state 0 = m_PS ground
#     jk     : binsize 10 (autocorrelation)
#   Overlays 2*m_PS(fit) against the two-meson m_1 = 0.6145(58) from [0,3,6] 2-op.
#   Reuses hankel_off / staged_project / rebased_effmass_fixed from hankel_rebase_scan_claude.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs

T0 = 3
REBT = 4
OFF = [0, 3, 6]
BINSIZE = 10
M1 = 0.6145        # two-meson state 1 from [0,3,6] 2-op
M1E = 0.0058
SCRATCH = hs.SCRATCH


def build_mps_store():
    tag = dc.ENS.split("nu0")[0]
    ncfg = len(dc.KS)
    cache = "%s/mps_store_%s_%d_claude.pkl" % (SCRATCH, tag.replace(".", "p"), ncfg)
    if os.path.exists(cache):
        with open(cache, "rb") as f:
            allM, twin = pickle.load(f)
        print("# loaded cached m_PS store  ncfg=%d twin=%d" % (ncfg, twin))
        return allM, twin, ncfg, tag
    print("# building single-sigma M(t) store  ncfg=%d ..." % ncfg)
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allM = []
    twin = None
    for j, k in enumerate(dc.KS):
        V, tau, taugw, tsrc0, tw = dc.load_peram(k)
        twin = tw
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(tw)]
        M = np.zeros(tw)
        for dt in range(tw):
            ns = tw - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                acc += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            M[dt] = acc / ns
        allM.append(M)
        if (j + 1) % 100 == 0:
            print("#   ... %d/%d" % (j + 1, ncfg))
    allM = np.array(allM)
    with open(cache, "wb") as f:
        pickle.dump((allM, twin), f)
    print("# cached m_PS store -> %s" % cache)
    return allM, twin, ncfg, tag


def hankel_effmass(Mvec, twin, Vfix):
    Cts = Mvec.reshape(twin, 1, 1)
    Big = hs.hankel_off(Cts, OFF)
    if Vfix is None:
        V = hs.staged_project(Big, [(REBT, 2)], T0)
        return hs.rebased_effmass_fixed(Big, V, T0), V
    return hs.rebased_effmass_fixed(Big, Vfix, T0), Vfix


allM, twin, ncfg, tag = build_mps_store()
nbin = ncfg // BINSIZE
blocks = np.array([allM[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nbin)])

em_c, Vfix = hankel_effmass(allM.mean(0), twin, None)
ems = []
for i in range(nbin):
    Mi = np.delete(blocks, i, 0).mean(0)
    em_i, _ = hankel_effmass(Mi, twin, Vfix)
    ems.append(em_i)
ems = np.array(ems)
em_err = np.sqrt((nbin - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
tmax = em_c.shape[0]

print("\n#  t |   mPS0(err)       mPS1(err)")
for t in range(T0, min(tmax, 20)):
    print("#  %2d | %7.4f(%.4f)  %7.4f(%.4f)" % (t, em_c[t, 0], em_err[t, 0], em_c[t, 1], em_err[t, 1]))


def corr_const_fit(lev, tlo, thi):
    ts = np.arange(tlo, thi + 1)
    y = em_c[ts, lev]
    Y = ems[:, ts, lev]
    Ym = Y.mean(0)
    Cov = (nbin - 1.0) / nbin * np.einsum("ip,iq->pq", Y - Ym, Y - Ym)
    Cinv = np.linalg.inv(Cov)
    one = np.ones(len(ts))
    denom = one @ Cinv @ one
    mhat = (one @ Cinv @ y) / denom
    err = np.sqrt(1.0 / denom)
    r = y - mhat
    return mhat, err, r @ Cinv @ r, len(ts) - 1


# choose a plateau window for the ground; report a couple
for (a, b) in [(6, 9), (8, 11), (10, 13)]:
    m, e, c2, dof = corr_const_fit(0, a, b)
    print("# m_PS fit [%2d,%2d]: %.4f(%d)  chi2/dof=%.2f/%d  ->  2m_PS=%.4f(%d)"
          % (a, b, m, round(e * 1e4), c2, dof, 2 * m, round(2 * e * 1e4)))

mfit, efit, _, _ = corr_const_fit(0, 8, 11)
two = 2 * mfit
twoe = 2 * efit
print("\n# 2m_PS(clean) = %.4f(%d)   vs   m_1(two-meson) = %.4f(%d)"
      % (two, round(twoe * 1e4), M1, round(M1E * 1e4)))
print("# Delta = m_1 - 2m_PS = %.4f  (%.1f sigma)"
      % (M1 - two, (M1 - two) / np.sqrt(twoe ** 2 + M1E ** 2)))

# ---- plot ----
ts = np.arange(tmax)
fig, ax = plt.subplots(figsize=(8.6, 5.6))
g0 = np.isfinite(em_c[:, 0]) & (em_err[:, 0] < 0.2)
g1 = np.isfinite(em_c[:, 1]) & (em_err[:, 1] < 0.2)
ax.errorbar(ts[g0], em_c[g0, 0], yerr=em_err[g0, 0], color="tab:green", marker="o", ms=5, lw=1.1,
            capsize=2.5, label=r"$m_{PS}$ ground (state 0)")
ax.errorbar(ts[g1], em_c[g1, 1], yerr=em_err[g1, 1], color="tab:purple", marker="v", ms=5, lw=1.1,
            capsize=2.5, alpha=0.7, label=r"meson excited (state 1)")
ax.fill_between([T0, min(tmax, 20)], two - twoe, two + twoe, color="tab:blue", alpha=0.2)
ax.plot([T0, min(tmax, 20)], [two, two], color="tab:blue", lw=1.6, label=r"$2m_{PS}$ (clean fit)")
ax.fill_between([T0, min(tmax, 20)], M1 - M1E, M1 + M1E, color="tab:red", alpha=0.2)
ax.plot([T0, min(tmax, 20)], [M1, M1], color="tab:red", lw=1.6, ls="--",
        label=r"$m_1$ two-meson $=0.6145(58)$")
ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.4)
ax.set_ylim(0.2, 0.9)
ax.set_xlim(T0, min(tmax, 20))
ax.set_xlabel(r"$t$")
ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
ax.set_title("single-$\\sigma$ $m_{PS}$ via Hankel [0,3,6] 2-op T0=3 reb 2@4  vs  two-meson $m_1$", fontsize=11)
ax.legend(fontsize=9, loc="upper right")
ax.grid(alpha=0.3)
fig.tight_layout()
out = "figs/hankel_mPS_dt036_vs_m1_claude.png"
os.makedirs("figs", exist_ok=True)
fig.savefig(out, dpi=130)
plt.close(fig)
print("# -> %s" % out)
