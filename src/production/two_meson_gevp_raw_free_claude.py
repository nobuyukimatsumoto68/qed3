#!/usr/bin/env python3
# two_meson_gevp_raw_free_claude.py  [FREE -- {1, sigma_00, sigma_00^2} wall (l=0) GEVP]
# Run:  ENS=free LREF=2 NVDIR=distill_Nv84 CONTACT=0 python3 two_meson_gevp_raw_free_claude.py
#
# UNIMPROVED overlap propagator: sigma = psibar D_ov^{-1} psi with NO contact subtraction (raw tau).
# CONTACT env sets the equal-time subtraction: 0.0 = unimproved (this run), 0.5 = improved (tilde_tau).
# {1, sigma_00, sigma_00^2} 3x3 GEVP, identity carries the vacuum (its off-diagonals are the measured
# equal-time one-points <sigma>, <sigma^2>).  All blocks via the general mode-space Wick engine.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from itertools import permutations
import distill_contract_claude as dc

T0 = 3
CONTACT = float(os.environ.get("CONTACT", "0.0"))     # 0 = unimproved (raw), 0.5 = improved


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


def solve_gevp(Cts, T0, tol=1e-10):
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
    print("# ENS=%s  FREE  {1, sigma_00, sigma_00^2} wall GEVP  t0=%d  L=%d  CONTACT=%.2f (%s)"
          % (tag, T0, dc.L, CONTACT, "unimproved/raw" if CONTACT == 0.0 else "improved"))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00

    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - CONTACT * Iv                # equal-time subtraction (0 = raw)

    P00 = []
    for a in range(twin):
        Vt = V[tsrc0 + a].T
        P00.append(Vt.conj().T @ (w00[:, None] * Vt))

    # equal-time one-points (identity off-diagonals)
    o1 = np.mean([wick([(P00[a], a)], tt).real for a in range(twin)])                 # <sigma>
    o2 = np.mean([wick([(P00[a], a), (P00[a], a)], tt).real for a in range(twin)])    # <sigma^2>
    print("# one-points: <sigma> = %.6e   <sigma^2> = %.6e" % (o1, o2))

    C11 = np.zeros(twin)
    C12 = np.zeros(twin)
    C22 = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a11 = a12 = a22 = 0.0
        for s in range(ns):
            t = s + dt
            a11 += wick([(P00[t], t), (P00[s], s)], tt).real
            a12 += wick([(P00[t], t), (P00[s], s), (P00[s], s)], tt).real
            a22 += wick([(P00[t], t), (P00[t], t), (P00[s], s), (P00[s], s)], tt).real
        C11[dt], C12[dt], C22[dt] = a11 / ns, a12 / ns, a22 / ns
    print("# large-t: C11=%.4e C12=%.4e C22=%.4e  (vac: C22->o2^2=%.4e, C11->o1^2=%.4e)"
          % (C11[twin - 3], C12[twin - 3], C22[twin - 3], o2 ** 2, o1 ** 2))

    Cts = np.zeros((twin, 3, 3))
    for dt in range(twin):
        Cts[dt] = np.array([[1.0, o1, o2],
                            [o1, C11[dt], C12[dt]],
                            [o2, C12[dt], C22[dt]]])

    lam, nlev = solve_gevp(Cts, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    msig = {1: 0.378, 2: 0.393}.get(dc.L, None)
    print("# levels kept = %d ; m_sigma(L%d) ~ %s ; 2 m_sigma ~ %s"
          % (nlev, dc.L, msig, None if msig is None else 2 * msig))
    print("\n  t   " + "  ".join("m%d" % i for i in range(nlev)))
    for t in range(T0 + 1, twin - 2):
        row = "  ".join("%7.4f" % em[t, i] if np.isfinite(em[t, i]) else "   nan " for i in range(nlev))
        print("  %2d  %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % len(cols)], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, msig or 0.0, 2 * (msig or 0.0)]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE {$1,\sigma,\sigma^2$} wall GEVP (L=%d, CONTACT=%.2f)" % (dc.L, CONTACT))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_raw_free_L%d_c%.1f_claude.png" % (dc.L, CONTACT)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
