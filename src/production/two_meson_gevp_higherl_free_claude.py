#!/usr/bin/env python3
# two_meson_gevp_higherl_free_claude.py  [FREE -- {1, sigma_00, sigma_00^2, sigma_{J6,J10,J12}} GEVP]
# Run:  ENS=free LREF=2 NVDIR=distill_Nv84 python3 two_meson_gevp_higherl_free_claude.py
#
# Adds icosahedral-invariant HIGHER-l single-meson operators sigma_{J_l} = sum_x w_x J_l(x) sigma(x)
# for l = 6, 10, 12 (the lowest orbital harmonics that subduce to the trivial A irrep on the icos
# lattice; icos_harmonics_claude.py).  J_l = continuum icosahedral harmonic (reusing the same Y_lm
# convention as s2.h: m = 0 mod 5, theta=acos z, phi=atan2 y,x), sampled at the mesh sites.
# NOTE: the mesh supports only (#orbits) invariant operators -- L2->1, L3->2, L4->3 -- so at coarse L
# the higher-l set is rank-deficient; the GEVP solver prunes the null space of C(t0) automatically.
# All blocks via the general mode-space Wick engine (uniform normalization).  Identity carries the vacuum.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from itertools import permutations
import distill_contract_claude as dc
import icos_harmonics_claude as ic

T0 = 3
LVALS = [6, 10, 12]


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
    # generalized eigenvalues of (C(t), C(t0)) with C(t0) null space pruned (whitening on its range)
    twin = Cts.shape[0]
    C0 = 0.5 * (Cts[T0] + Cts[T0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > tol * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])          # whitening: Uk^T C0 Uk = I
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
    print("# ENS=%s  FREE  {1, sigma_00, sigma_00^2, sigma_{J6,J10,J12}} GEVP  t0=%d  L=%d"
          % (tag, T0, dc.L))
    dual = dc.dual_areas_from_mesh()
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    sites = sites / np.linalg.norm(sites, axis=1, keepdims=True)
    verts = dc.load_vec3(dc.GEOM + "pts_n1.dat")
    verts = verts / np.linalg.norm(verts, axis=1, keepdims=True)
    rots = ic.icosahedral_rotations(verts)

    # vertex spatial weights: index 0 = W00 (l=0), then J_l for l in LVALS
    wvecs = [np.repeat(dual, dc.NS) * dc.Y00]
    labels = ["s00"]
    for l in LVALS:
        Jl = ic.Jl_on_sites(sites, l, rots)
        wvecs.append(np.repeat(dual * Jl, dc.NS))
        labels.append("J%d" % l)

    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv

    # vertices P[k][a] (Nv x Nv) for each spatial weight k, each timeslice a
    P = []
    for w in wvecs:
        Pk = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Pk.append(Vt.conj().T @ (w[:, None] * Vt))
        P.append(Pk)

    # operator descriptors: list of vertex-indices (legs).  identity handled separately.
    #   sigma_00 = [0] ; sigma_00^2 = [0,0] ; sigma_{J_l} = [k]
    ops = [("s00", [0]), ("s00^2", [0, 0])]
    for k in range(1, len(wvecs)):
        ops.append((labels[k], [k]))
    nop = len(ops)                                  # non-identity operators
    print("# non-identity operators: %s" % ", ".join(o[0] for o in ops))

    # equal-time one-points <op> (identity off-diagonals): single-leg ops -> ~0; sigma_00^2 -> o2
    def onepoint(desc):
        val = 0.0
        for a in range(twin):
            legs = [(P[k][a], a) for k in desc]
            val += wick(legs, tt).real
        return val / twin
    ovec = np.array([onepoint(d) for (_, d) in ops])
    print("# one-points <op> = %s" % np.array2string(ovec, precision=4))

    # correlator blocks C[dt][i,j] for non-identity ops (sink at t, source at s)
    C = np.zeros((twin, nop, nop))
    for dt in range(twin):
        ns = twin - dt
        acc = np.zeros((nop, nop))
        for s in range(ns):
            t = s + dt
            for i in range(nop):
                for j in range(i, nop):
                    legs = [(P[k][t], t) for k in ops[i][1]] + [(P[k][s], s) for k in ops[j][1]]
                    v = wick(legs, tt).real
                    acc[i, j] += v
                    if j != i:
                        acc[j, i] += v
        C[dt] = acc / ns

    # assemble full matrices with the identity as operator 0
    n = nop + 1
    Cts = np.zeros((twin, n, n))
    for dt in range(twin):
        Cts[dt, 0, 0] = 1.0
        Cts[dt, 0, 1:] = ovec
        Cts[dt, 1:, 0] = ovec
        Cts[dt, 1:, 1:] = C[dt]

    lam, nlev = solve_gevp(Cts, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    mline = {1: 0.393, 2: 0.786}.get(dc.L, None)
    print("# GEVP levels kept (rank of C(t0)) = %d" % nlev)
    print("\n  t   " + "  ".join("m%d" % i for i in range(nlev)))
    for t in range(T0 + 1, twin - 2):
        row = "  ".join("%6.3f" % em[t, i] if np.isfinite(em[t, i]) else "  nan " for i in range(nlev))
        print("  %2d  %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue", "tab:green", "tab:purple", "tab:orange"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % len(cols)], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, 0.393, 0.786]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.6)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE {$1,\sigma_{00},\sigma_{00}^2,\sigma_{J_\ell}$} GEVP (L=%d, %d levels)" % (dc.L, nlev))
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_higherl_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
