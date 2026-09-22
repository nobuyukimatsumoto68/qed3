#!/usr/bin/env python3
# two_meson_gevp_pp_free_claude.py  [FREE -- {1, sigma_00, sigma_00^2, O^(2)_1} 4x4 GEVP]
# Run:  ENS=free LREF=2 NVDIR=distill_Nv84 python3 two_meson_gevp_pp_free_claude.py
#
# Adds the relative-l=1 back-to-back two-meson operator O^(2)_1 = sum_mu sigma_{1,mu}^2 (Cartesian
# real-vector L=0 singlet, sigma_pp_operator_note_claude.md) to the {1, sigma_00, sigma_00^2} basis.
# All 4-, 3-, 2-point blocks via a GENERAL mode-space Wick sum (analogue of the validated wick4 in
#   point_source_full_claude.py):  <prod_i psibar Phi_i psi> = sum_{pi in S_n} prod_cycles (-1) Tr[..],
#   G_ij = Phi_i tilde_tau(t_i,t_j),  tilde_tau(a,b) = tau(a,b) - (1/2) delta_ab I  (contact-subtracted).
# Identity operator carries the vacuum: its off-diagonals are the measured equal-time one-points
#   <sigma_00^2>, <O_1>.  Level 0 -> vacuum (0), level 1 -> m_sigma, and we test whether a two-meson
# level near 2 m_sigma separates from the excited-single band.
# Refs: Peardon 0905.2160 (distillation); Luscher-Wolff (1990), Blossier 0902.1265 (variational).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from itertools import permutations
from scipy.linalg import eig
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = 3
NORM1 = np.sqrt(3.0 / (4.0 * np.pi))              # l=1 real harmonic normalization (overall scale only)


def wick(legs, tt):
    # legs: list of (Phi, tidx).  tt: tilde_tau tensor (twin,twin,Nv,Nv), diagonal already contact-subtracted.
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
    de.CONTACT = 0.0
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE  {1, sigma_00, sigma_00^2, O^(2)_1} 4x4 GEVP  t0=%d  L=%d" % (tag, T0, dc.L))
    dual = dc.dual_areas_from_mesh()
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)          # (N_sites, 3) unit-sphere coords
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    w1 = [np.repeat(dual * sites[:, mu], dc.NS) * NORM1 for mu in range(3)]   # l=1 vertices, mu=x,y,z

    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv                          # tilde_tau on equal-time diagonal

    # vertices per timeslice
    P00 = []
    P1 = [[], [], []]
    for a in range(twin):
        Vt = V[tsrc0 + a].T                                      # (2Ns, Nv), column k = mode k
        Vd = Vt.conj().T
        P00.append(Vd @ (w00[:, None] * Vt))
        for mu in range(3):
            P1[mu].append(Vd @ (w1[mu][:, None] * Vt))

    # equal-time one-points (identity off-diagonals):  <sigma_a^2> = <sigma_a sigma_a>, tilde_tau
    o2 = np.mean([wick([(P00[a], a), (P00[a], a)], tt).real for a in range(twin)])
    oO = np.mean([sum(wick([(P1[mu][a], a), (P1[mu][a], a)], tt).real for mu in range(3))
                  for a in range(twin)])
    print("# one-points: <sigma_00^2> = %.6e   <O_1> = %.6e" % (o2, oO))

    C11 = np.zeros(twin)
    C12 = np.zeros(twin)
    C22 = np.zeros(twin)
    C1O = np.zeros(twin)
    C2O = np.zeros(twin)
    COO = np.zeros(twin)
    # ALL blocks via the one Wick engine -> uniform normalization (sink at t, source at s).
    # C22 here = wick4 = sum_i W10_i diags_i (the TRUE correlator); the old closed form carried a
    # spurious extra x2 (sum W10 = 24 = |S_4| already gives the full contraction), which would be
    # inconsistent with C11, C12 and the l=1 blocks.  We keep the closed forms as cross-checks only.
    for dt in range(twin):
        ns = twin - dt
        a11 = a12 = a22 = a1O = a2O = aOO = 0.0
        for s in range(ns):
            t = s + dt
            a11 += wick([(P00[t], t), (P00[s], s)], tt).real
            a12 += wick([(P00[t], t), (P00[s], s), (P00[s], s)], tt).real
            a22 += wick([(P00[t], t), (P00[t], t), (P00[s], s), (P00[s], s)], tt).real
            for mu in range(3):
                a1O += wick([(P00[t], t), (P1[mu][s], s), (P1[mu][s], s)], tt).real
                a2O += wick([(P00[t], t), (P00[t], t), (P1[mu][s], s), (P1[mu][s], s)], tt).real
                for nu in range(3):
                    aOO += wick([(P1[mu][t], t), (P1[mu][t], t),
                                 (P1[nu][s], s), (P1[nu][s], s)], tt).real
        C11[dt], C12[dt], C22[dt] = a11 / ns, a12 / ns, a22 / ns
        C1O[dt], C2O[dt], COO[dt] = a1O / ns, a2O / ns, aOO / ns

    # cross-check vs the old closed forms at one dt
    dtc = min(6, twin - 1)
    s = 0
    tc = s + dtc
    c11_cf = (-np.trace(P00[s] @ tau[s, tc] @ P00[tc] @ tau[tc, s])).real
    m = P00[s] @ tau[s, tc] @ P00[tc] @ tau[tc, s]
    c12_cf = (-2.0 * np.trace(P00[s] @ tt[s, s] @ m)).real
    c22_cf = 2.0 * (dc.W10 * de.diags_pair(P00, tt, s, tc)).sum().real
    print("# CHECK dt=%d:  C11 wick=%.4e cf=%.4e | C12 wick=%.4e cf=%.4e | C22 wick=%.4e cf(old x2)=%.4e"
          % (dtc, C11[dtc], c11_cf, C12[dtc], c12_cf, C22[dtc], c22_cf))
    print("# large-t: sqrt(C22)=%.4e vs <sigma^2>=%.4e ; sqrt(COO)=%.4e vs <O_1>=%.4e"
          % (np.sqrt(abs(C22[twin - 3])), o2, np.sqrt(abs(COO[twin - 3])), oO))

    def Cmat(t):
        return np.array([[1.0, 0.0, o2, oO],
                         [0.0, C11[t], C12[t], C1O[t]],
                         [o2, C12[t], C22[t], C2O[t]],
                         [oO, C1O[t], C2O[t], COO[t]]])

    nlev = 4
    lam = np.full((twin, nlev), np.nan)
    C0 = Cmat(T0)
    for t in range(twin):
        try:
            lam[t] = np.sort(eig(Cmat(t), C0)[0].real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])

    print("\n  t    m0       m1       m2       m3     [0 ; 0.393 ; 0.786 (L2)]")
    for t in range(T0 + 1, twin - 2):
        row = "  ".join("%7.4f" % em[t, i] if np.isfinite(em[t, i]) else "   nan " for i in range(nlev))
        print("  %2d   %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue", "tab:green"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, 0.393, 0.786]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE {$1,\sigma_{00},\sigma_{00}^2,O^{(2)}_1$} 4x4 GEVP (L=%d)" % dc.L)
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_pp_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
