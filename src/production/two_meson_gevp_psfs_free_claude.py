#!/usr/bin/env python3
# two_meson_gevp_psfs_free_claude.py  [FREE -- {1, sigma_PS, sigma_PS^2, sigma_FS} 4x4 GEVP]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_gevp_psfs_free_claude.py
#
# Add the FURNISHED scalar sigma_FS (single-meson leg = -taugw = -(1-D_ov^dag)D^{-1}) to the PS basis.
# Motivation: <sigma_FS sigma_PS^2> is a nonzero, smooth 2->1 amplitude (two_meson_ps2_to_fs_free), while
# <sigma_PS sigma_PS^2> ~ 0.  Does the FS mixing let the GEVP resolve a two-meson level?
# sigma_FS DEFINITION here = pure furnished operator (leg gw); NOT the codebase Vpp+Vff combination -- easy
# to switch if wanted.  Blocks by explicit traces; C22 (PS two-meson) via the Wick engine.  Contact tilde_tau
# on equal-time PS legs.  Identity carries the vacuum (measured one-points).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from itertools import permutations
import distill_contract_claude as dc

T0 = 3


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


def solve_gevp(Cts, T0, tol=1e-11):
    twin = Cts.shape[0]
    C0 = 0.5 * (Cts[T0] + Cts[T0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > tol * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    nlev = int(keep.sum())
    lam = np.full((twin, nlev), np.nan)
    for t in range(twin):
        M = Uk.T @ (0.5 * (Cts[t] + Cts[t].T)) @ Uk
        try:
            lam[t] = np.sort(np.linalg.eigvals(M).real)[::-1]
        except Exception:
            pass
    return lam, nlev


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  {1, sigma_PS, sigma_PS^2, sigma_FS} 4x4 GEVP  t0=%d  m_sigma~%.3f"
          % (tag, dc.L, T0, msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv                # tilde_tau
    gw = -taugw                                        # FS furnished leg
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    # one-points (identity off-diagonals), equal-time, mean over time
    oPS = np.mean([(-np.trace(Phi[a] @ tt[a, a])).real for a in range(twin)])
    oFS = np.mean([(-np.trace(Phi[a] @ gw[a, a])).real for a in range(twin)])
    oPS2 = np.mean([wick([(Phi[a], a), (Phi[a], a)], tt).real for a in range(twin)])
    print("# one-points: <sigma_PS>=%.4e  <sigma_FS>=%.4e  <sigma_PS^2>=%.4e" % (oPS, oFS, oPS2))

    # correlator blocks (sink t, source s = t - dt), translation-averaged
    P = np.zeros(twin)          # <PS PS>
    Qff = np.zeros(twin)        # <FS FS>
    Rpf = np.zeros(twin)        # <PS FS> (symmetrized below)
    C22 = np.zeros(twin)        # <PS^2 PS^2>
    triPP = np.zeros(twin)      # <PS PS^2>
    triFS = np.zeros(twin)      # <FS PS^2>
    for dt in range(twin):
        ns = twin - dt
        aP = aQ = aR = aC = aTp = aTf = 0.0
        for s in range(ns):
            t = s + dt
            aP += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            aQ += (-np.trace(Phi[t] @ gw[t, s] @ Phi[s] @ gw[s, t])).real
            aR += (-0.5 * (np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ gw[s, t])
                           + np.trace(Phi[t] @ gw[t, s] @ Phi[s] @ tau[s, t]))).real
            aC += wick([(Phi[t], t), (Phi[t], t), (Phi[s], s), (Phi[s], s)], tt).real
            # triangles, symmetrized over which time holds sigma^2
            aTp += (-np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s])
                    - np.trace(Phi[t] @ tt[t, t] @ Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            aTf += (-np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ gw[t, s])
                    - np.trace(Phi[t] @ tt[t, t] @ Phi[t] @ tau[t, s] @ Phi[s] @ gw[s, t])).real
        P[dt], Qff[dt], Rpf[dt] = aP / ns, aQ / ns, aR / ns
        C22[dt], triPP[dt], triFS[dt] = aC / ns, aTp / ns, aTf / ns

    # PLAIN {1, sigma_FS, sigma_PS^2} GEVP -- identity carries the vacuum (drop single sigma_PS).
    # Qff and triFS were computed CONNECTED only; add their DISCONNECTED vacuum pieces so every block
    # C_XY -> <X><Y> at large t (C22 already includes its disconnected F=D'D' term).  order [1, FS, PS^2]
    Cts = np.zeros((twin, 3, 3))
    for dt in range(twin):
        Cts[dt] = np.array([[1.0,   oFS,              oPS2],
                            [oFS,   Qff[dt] + oFS**2,  triFS[dt] + oFS * oPS2],
                            [oPS2,  triFS[dt] + oFS * oPS2, C22[dt]]])

    lam, nlev = solve_gevp(Cts, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# PLAIN {1, sigma_FS, sigma_PS^2} GEVP ; levels kept = %d ; m_sigma=%.3f 2m_sigma=%.3f"
          % (nlev, msig, 2 * msig))

    print("\n  t  | {1, FS, PS^2}: m0      m1      m2")
    for t in range(T0 + 1, twin - 2):
        r3 = "  ".join("%7.4f" % em[t, i] if i < nlev and np.isfinite(em[t, i]) else "   --- " for i in range(3))
        print("  %2d | %s" % (t, r3))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue", "tab:green"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % 4], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, msig, 2 * msig]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d plain {$1,\sigma_{FS},\sigma_{PS}^2$} GEVP" % dc.L)
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_psfs_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
