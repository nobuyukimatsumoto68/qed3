#!/usr/bin/env python3
# two_meson_mix_diagrams_free_claude.py  [FREE -- 1<->2 (sigma <-> sigma^2) triangle, diagram-by-diagram]
# Run:  ENS=free LREF=2 NVDIR=distill_Nv84 python3 two_meson_mix_diagrams_free_claude.py
#
# The 1<->2 mixing correlator (NM contraction_v3.pdf) is  <sigma^2 sigma> = 2C' + 2G' + I' + J', where
#   C'  = -Tr[Phi(0) tilde_tau(0,0) Phi(0) tau(0,t) Phi(t) tau(t,0)]      (the connected TRIANGLE)
#   G'  = D_S . C_S ,  I' = D'_S . D_S(t) ,  J' = D_S^2 . D_S(t)          (tadpoles: standalone D_S -> 0 under tilde_tau)
# We compute BOTH directions (2->1: sigma^2 at source; 1->2: sigma^2 at sink), show the triangle
# correlator + effmass and its smoothness/error floor, and confirm the tadpole remainder ~ 0.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from itertools import permutations
import distill_contract_claude as dc


def wick(legs, tt):
    n = len(legs)
    G = [[legs[i][0] @ tt[legs[i][1], legs[j][1]] for j in range(n)] for i in range(n)]
    total = 0.0 + 0.0j
    for perm in permutations(range(n)):
        visited = [False] * n
        val = 1.0 + 0.0j
        for start in range(n):
            if visited[start]:
                continue
            i = start
            prod = None
            while not visited[i]:
                visited[i] = True
                blk = G[i][perm[i]]
                prod = blk if prod is None else prod @ blk
                i = perm[i]
            val *= (-1.0) * np.trace(prod)
        total += val
    return total


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  1<->2 (sigma<->sigma^2) triangle diagrams  m_sigma~%.3f" % (tag, dc.L, msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv                        # tilde_tau
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    # per dt, translation-averaged.  Direction 2->1: sigma^2 at source s, sigma at sink t=s+dt.
    #               Direction 1->2: sigma at source s, sigma^2 at sink t.
    tri21 = np.zeros(twin); full21 = np.zeros(twin)
    tri12 = np.zeros(twin); full12 = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a21 = f21 = a12 = f12 = 0.0
        for s in range(ns):
            t = s + dt
            # 2->1 triangle:  -2 Tr[Phi_s tilde_tau(s,s) Phi_s tau(s,t) Phi_t tau(t,s)]
            m = Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s]
            a21 += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ m)).real
            f21 += wick([(Phi[t], t), (Phi[s], s), (Phi[s], s)], tt).real     # full 2->1 (sigma^2 at s)
            # 1->2 triangle:  -2 Tr[Phi_t tilde_tau(t,t) Phi_t tau(t,s) Phi_s tau(s,t)]
            m2 = Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t]
            a12 += (-2.0 * np.trace(Phi[t] @ tt[t, t] @ m2)).real
            f12 += wick([(Phi[t], t), (Phi[t], t), (Phi[s], s)], tt).real     # full 1->2 (sigma^2 at t)
        tri21[dt], full21[dt] = a21 / ns, f21 / ns
        tri12[dt], full12[dt] = a12 / ns, f12 / ns
    tad21 = full21 - tri21
    tad12 = full12 - tri12

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))

    print("\n  dt |  tri(2->1)    full(2->1)   tad(2->1)  |  tri(1->2)    tad(1->2)  | m_eff tri21  tri12")
    e21 = eff(tri21); e12 = eff(tri12)
    for dt in range(1, min(24, twin - 1)):
        print("  %2d | %+.3e  %+.3e  %+.2e | %+.3e  %+.2e |  %6.3f     %6.3f"
              % (dt, tri21[dt], full21[dt], tad21[dt], tri12[dt], tad12[dt], e21[dt], e12[dt]))
    print("# max|tadpole remainder| : 2->1 = %.2e , 1->2 = %.2e  (should be ~0 under tilde_tau)"
          % (np.abs(tad21).max(), np.abs(tad12).max()))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)

    # panel 1: |triangle| correlator (log), both directions + tadpole remainder
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    ax.semilogy(dts, np.abs(tri21[1:]), color="tab:red", marker="o", ms=3, lw=1, label=r"triangle $2C'$ (2$\to$1)")
    ax.semilogy(dts, np.abs(tri12[1:]), color="tab:blue", marker="s", ms=3, lw=1, label=r"triangle $2C'$ (1$\to$2)")
    ax.semilogy(dts, np.abs(tad21[1:]), color="tab:gray", marker="^", ms=3, lw=1, label="tadpole remainder (2$\\to$1)")
    ax.set_xlabel("dt"); ax.set_ylabel(r"$|C(dt)|$")
    ax.set_title(r"FREE L=%d  1$\leftrightarrow$2 triangle: correlator" % dc.L)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out1 = "figs/two_meson_mix_corr_free_L%d_claude.png" % dc.L
    os.makedirs("figs", exist_ok=True)
    fig.savefig(out1, dpi=130)
    plt.close(fig)

    # panel 2: triangle effmass, both directions
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(tri21, "tab:red", "o", r"$2C'$ (2$\to$1)"), (tri12, "tab:blue", "s", r"$2C'$ (1$\to$2)")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=3, lw=1, label=lab)
    for y in [msig, 2 * msig]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("dt"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  1$\leftrightarrow$2 triangle: effmass (dashed $m_\sigma$, $2m_\sigma$)" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out2 = "figs/two_meson_mix_effmass_free_L%d_claude.png" % dc.L
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print("\n# -> %s\n# -> %s" % (out1, out2))


if __name__ == "__main__":
    main()
