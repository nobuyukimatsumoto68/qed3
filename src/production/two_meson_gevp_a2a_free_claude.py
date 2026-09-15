#!/usr/bin/env python3
# two_meson_gevp_a2a_free_claude.py  [FREE L1 -- {1, sigma_00, sigma_00^2} wall GEVP: ALL-TO-ALL vs distillation]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_gevp_a2a_free_claude.py
#
# At L1 the distillation basis is COMPLETE (Nv = 2*Ns = 24), so we can reconstruct the FULL position-space
# propagator G(x,t;y,t') = V(t)^T tau(t,t') V(t')^*  (24x24 per time-pair, "all-to-all", same cost) and
# contract the wall operators DIRECTLY in position space -- bypassing the mode-space (distillation) traces.
# We run BOTH and print them side by side; if distillation is exact they must agree.  W = diag_x(w_x Y00).
# Contact = -1/2 I on the equal-time diagonal (spatial delta in position space; improved propagator).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from itertools import permutations
import distill_contract_claude as dc

T0 = 3
CONTACT = float(os.environ.get("CONTACT", "0.5"))


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


def solve_gevp(Cts, T0, tol=1e-12):
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


def blocks(P, tt, twin):
    # {sigma, sigma^2} correlators + one-points, generic in the (vertex Phi list P, propagator tt)
    o1 = np.mean([wick([(P[a], a)], tt).real for a in range(twin)])
    o2 = np.mean([wick([(P[a], a), (P[a], a)], tt).real for a in range(twin)])
    C11 = np.zeros(twin)
    C12 = np.zeros(twin)
    C22 = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a11 = a12 = a22 = 0.0
        for s in range(ns):
            t = s + dt
            a11 += wick([(P[t], t), (P[s], s)], tt).real
            a12 += wick([(P[t], t), (P[s], s), (P[s], s)], tt).real
            a22 += wick([(P[t], t), (P[t], t), (P[s], s), (P[s], s)], tt).real
        C11[dt], C12[dt], C22[dt] = a11 / ns, a12 / ns, a22 / ns
    Cts = np.zeros((twin, 3, 3))
    for dt in range(twin):
        Cts[dt] = np.array([[1.0, o1, o2], [o1, C11[dt], C12[dt]], [o2, C12[dt], C22[dt]]])
    return Cts, o1, o2


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE L=%d  {1,sigma,sigma^2} wall GEVP  ALL-TO-ALL vs distillation  CONTACT=%.2f"
          % (tag, dc.L, CONTACT))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00                       # length 2Ns
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = tau.shape[-1]
    twoNs = V.shape[2]
    print("# Nv=%d  2Ns=%d  (complete basis: Nv==2Ns is %s)" % (Nv, twoNs, Nv == twoNs))
    Imode = np.eye(Nv)
    Ipos = np.eye(twoNs)

    # ---- distillation (mode-space) ----
    ttm = tau.copy()
    for a in range(twin):
        ttm[a, a] = tau[a, a] - CONTACT * Imode
    Pmode = []
    for a in range(twin):
        Vt = V[tsrc0 + a].T                                    # (2Ns, Nv)
        Pmode.append(Vt.conj().T @ (w00[:, None] * Vt))        # (Nv, Nv)
    Cm, o1m, o2m = blocks(Pmode, ttm, twin)

    # ---- all-to-all (position-space) ----  G(a,b) = V(a)^T tau(a,b) V(b)^*   (2Ns x 2Ns)
    Gpos = np.zeros((twin, twin, twoNs, twoNs), complex)
    for a in range(twin):
        Va = V[tsrc0 + a]                                       # (Nv, 2Ns)
        for b in range(twin):
            Vb = V[tsrc0 + b]
            Gpos[a, b] = Va.T @ tau[a, b] @ Vb.conj()
    for a in range(twin):
        Gpos[a, a] = Gpos[a, a] - CONTACT * Ipos               # spatial-delta contact
    Wpos = np.diag(w00)                                         # (2Ns, 2Ns) diagonal vertex
    Ppos = [Wpos for _ in range(twin)]
    Cp, o1p, o2p = blocks(Ppos, Gpos, twin)

    print("# one-points:  distill <sig>=%.6e <sig^2>=%.6e |  a2a <sig>=%.6e <sig^2>=%.6e"
          % (o1m, o2m, o1p, o2p))
    print("# max|C_distill - C_a2a| over all t = %.3e" % np.abs(Cm - Cp).max())

    lam_m, nlm = solve_gevp(Cm, T0)
    lam_p, nlp = solve_gevp(Cp, T0)
    with np.errstate(all="ignore"):
        em_m = np.log(lam_m[:-1] / lam_m[1:])
        em_p = np.log(lam_p[:-1] / lam_p[1:])
    msig = {1: 0.378, 2: 0.393}.get(dc.L, None)
    print("# levels kept: distill=%d  a2a=%d ; m_sigma(L%d)~%s" % (nlm, nlp, dc.L, msig))
    print("\n  t  | distill m0    m1     m2   | all-to-all m0   m1     m2")
    for t in range(T0 + 1, twin - 2):
        rm = "  ".join("%6.3f" % em_m[t, i] if i < nlm and np.isfinite(em_m[t, i]) else "  --- " for i in range(3))
        rp = "  ".join("%6.3f" % em_p[t, i] if i < nlp and np.isfinite(em_p[t, i]) else "  --- " for i in range(3))
        print("  %2d |  %s  |  %s" % (t, rm, rp))


if __name__ == "__main__":
    main()
