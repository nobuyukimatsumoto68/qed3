#!/usr/bin/env python3
# two_meson_prolif_compare_free_claude.py  -- does more block-Hankel proliferation do anything?
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_prolif_compare_free_claude.py
# NOVAC, t0=4, rebase into 2 states.  Overlays K time-shifts of {sigma^2, O_A}:
#   K=1 non-prolif (2 ops) ; K=2 shifts {0,1} (4 ops) ; K=3 shifts {0,1,2} (6 ops).
#   inflated correlator block (a,b) = c(t+a+b), a,b in 0..K-1 -> (2K)x(2K) block-Hankel.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = 4


def gevp_rebase(Cinf, tmax, nkeep, t0):
    C0 = 0.5 * (Cinf[t0] + Cinf[t0].T)
    wv, Uv = np.linalg.eigh(C0)
    order = np.argsort(wv)[::-1][:nkeep]
    Uk = Uv[:, order] / np.sqrt(wv[order])
    lam = np.full((Cinf.shape[0], nkeep), np.nan)
    for t in range(tmax):
        M = Uk.T @ (0.5 * (Cinf[t] + Cinf[t].T)) @ Uk
        try:
            lam[t] = np.sort(np.linalg.eigvals(M).real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    return em, lam


def inflate(c, K, twin):
    tmax = twin - 2 * (K - 1)
    Cinf = np.zeros((twin, 2 * K, 2 * K))
    for t in range(tmax):
        for a in range(K):
            for b in range(K):
                Cinf[t, 2 * a:2 * a + 2, 2 * b:2 * b + 2] = c[t + a + b]
    return Cinf, tmax


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    oo = np.outer([o2, oA], [o2, oA])

    C22 = np.zeros(twin); C2A = np.zeros(twin); CAA = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a22 = a2A = aAA = 0.0
        for s in range(ns):
            t = s + dt
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum().real
            a2A += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
        C22[dt], C2A[dt], CAA[dt] = a22 / ns, a2A / ns, aAA / ns
    c = np.zeros((twin, 2, 2))
    for t in range(twin):
        c[t] = np.array([[C22[t], C2A[t] + o2 * oA], [C2A[t] + o2 * oA, CAA[t] + oA ** 2]]) - oo   # connected

    res = {}
    for K in (1, 2, 3):
        Cinf, tmax = inflate(c, K, twin)
        em, lam = gevp_rebase(Cinf, tmax, 2, T0)
        res[K] = (em, lam, tmax)

    print("# NOVAC, t0=%d, rebase into 2.  two-meson (high) level for K=1,2,3 shifts:" % T0)
    print("\n  t  |  K=1     K=2     K=3   [2m_sig=0.756]")
    for t in range(T0 + 1, twin - 5):
        row = "  ".join("%6.3f" % res[K][0][t, 1] if t < res[K][2] - 1 and np.isfinite(res[K][0][t, 1]) else "  --- "
                        for K in (1, 2, 3))
        print("  %2d |  %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    styles = {1: ("tab:red", "o", "-"), 2: ("tab:green", "s", "--"), 3: ("tab:blue", "^", ":")}
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for K in (1, 2, 3):
        em, lam, tmax = res[K]
        col, mk, ls = styles[K]
        for i, lab in enumerate(["(2,2)", "two-meson"]):
            g = np.isfinite(em[:tmax - 1, i]) & (lam[:tmax - 1, i] * lam[1:tmax, i] > 0)
            ax.plot(ts[:tmax - 1][g], em[:tmax - 1][g, i], color=col, marker=mk, ms=4, lw=1.1, ls=ls,
                    label="K=%d %s" % (K, lab) if i == 1 else None)
    for y in [0.52, 2 * msig]:
        ax.axhline(y, color="k", ls=":", lw=0.8, alpha=0.5)
    ax.set_ylim(0.35, 0.95)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d NOVAC 2-state: proliferation K=1 (o), 2 (s), 3 (^)" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_prolif_compare_K_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
