#!/usr/bin/env python3
# sigma_ell_validate_propagator_claude.py
#   VALIDATE the continuum (Delta,ell,parity) classification against the FULL lattice propagator (tau).
#   Classification (sigma_pair_quantum_numbers, same/opposite-i3 split):
#     scalar   Gamma=1  (u*u+d*d): ell = |j1-j2|, +2, ...   ground: ell0=2E0, ell1=E0+E1, ell2=2E1/E0+E2
#     pseudo   Gamma=s3 (u*u-d*d): ell = |j1-j2|+-1, ...     ground: ell1(pseudo)=2E0
#   Here we compute the lattice correlator  C = -Tr[Phi(t) tau(t,s) Phi(s) tau(s,t)]  for each (Gamma,ell),
#   Phi = V^dag diag(Y_lm(x) (x) Gamma_spin) V, ell-multiplet M-summed, and read the effmass.
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 sigma_ell_validate_propagator_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import math
import numpy as np
import distill_contract_claude as dc


def real_harmonics(sites, l):
    x = sites[:, 0]
    y = sites[:, 1]
    z = sites[:, 2]
    if l == 0:
        return [np.full(len(x), 1.0 / math.sqrt(4.0 * math.pi))]
    if l == 1:
        c = math.sqrt(3.0 / (4.0 * math.pi))
        return [c * x, c * y, c * z]
    if l == 2:
        c1 = 0.25 * math.sqrt(5.0 / math.pi)
        c2 = 0.5 * math.sqrt(15.0 / math.pi)
        return [c2 * x * y, c2 * y * z, c1 * (3 * z * z - 1), c2 * x * z, 0.5 * c2 * (x * x - y * y)]


def main():
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    e0 = 0.5 * msig
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    dual = dc.dual_areas_from_mesh()
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    sites = sites / np.linalg.norm(sites, axis=1, keepdims=True)
    nsite = sites.shape[0]
    sign3 = np.tile([1.0, -1.0], nsite)                 # sigma3 on spin index j=NS*x+s
    Vt = [V[tsrc0 + a].T for a in range(twin)]

    def corr(wspin):
        Phi = [Vt[a].conj().T @ (wspin[:, None] * Vt[a]) for a in range(twin)]
        C = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                acc += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            C[dt] = acc / ns
        return C

    def channel(gamma, l):
        Yl = real_harmonics(sites, l)
        C = np.zeros(twin)
        for Y in Yl:
            w = np.repeat(dual * Y, dc.NS)
            if gamma == "s3":
                w = w * sign3
            C += corr(w)
        return C

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))

    channels = [("scalar", "id", 0, "2E0=%.3f" % (2 * e0)),
                ("scalar", "id", 1, "E0+E1~%.3f" % (e0 + 0.26)),
                ("scalar", "id", 2, "2E1/E0+E2~0.51"),
                ("pseudo", "s3", 0, "(none: no l=0)"),
                ("pseudo", "s3", 1, "2E0=%.3f" % (2 * e0)),
                ("pseudo", "s3", 2, "E0+E1~%.3f" % (e0 + 0.26))]
    Cs = {}
    print("# VALIDATION vs full propagator  FREE L=%d   2E0=%.3f  E0+E1~%.3f" % (dc.L, 2 * e0, e0 + 0.26))
    print("# channel        prediction         m_eff(dt=8)  (dt=15)  (dt=%d)" % (twin - 2))
    for (name, g, l, pred) in channels:
        C = channel(g, l)
        Cs[(name, l)] = C
        e = eff(C)
        d1, d2, d3 = 8, 15, twin - 2
        print("  %-7s l=%d  %-18s   %6.3f    %6.3f   %6.3f"
              % (name, l, pred, e[d1] if d1 < len(e) else np.nan,
                 e[d2] if d2 < len(e) else np.nan, e[d3] if d3 < len(e) else np.nan))

    print("\n# full effmass curves (scalar l0/l1/l2 ; pseudo l1):")
    print("#  dt | S_l0    S_l1    S_l2  | P_l1")
    for dt in range(2, twin - 1):
        row = (eff(Cs[("scalar", 0)])[dt], eff(Cs[("scalar", 1)])[dt],
               eff(Cs[("scalar", 2)])[dt], eff(Cs[("pseudo", 1)])[dt])
        print("  %2d | %6.3f  %6.3f  %6.3f | %6.3f" % (dt, *row))


if __name__ == "__main__":
    main()
