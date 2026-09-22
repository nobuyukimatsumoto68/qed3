#!/usr/bin/env python3
# e0e1_modes_continuum_claude.py
#   Pure-continuum enumeration: which local bilinear vertices Gamma = Y_lM(x) Sigma couple which
#   (m,n,i3) modes, and which have their LOWEST surviving fermion-antifermion pair at E0+E1
#   (lambda=1 x lambda=2) rather than the 2E0 ground.  See e0e1_continuum_impl_plan_claude.md.
#
#   Energies in continuum units E_k = lambda_k = n+|m|+1/2  (2E0<->2, E0+E1<->3, 2E1<->4).
#   Vertex matrix element  Gamma^lM_ab = int dOmega psi_a^dag Sigma psi_b Y_lM.
#   Channel strength (M-summed, rotation invariant):  P_l(a,b) = sum_M |Gamma^lM_ab|^2.
#   Free two-point:  C(t) = sum_ab P_l(a,b) exp(-(E_a+E_b) t)  -> effmass = smallest surviving E_a+E_b.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
from scipy.special import sph_harm_y
import free_wavefunctions_claude as fw

TOL = 1.0e-8                                  # |P_l| above this counts as "surviving"


def build_grid(N):
    gi = np.arange(N) + 0.5
    phi = np.mod(2 * np.pi * gi / ((1 + 5 ** 0.5) / 2), 2 * np.pi)
    ct = 1 - 2 * gi / N
    theta = np.arccos(ct)
    dV = 4 * np.pi / N
    return theta, phi, dV


def mode_label(m, n, i3):
    return "(m=%+.1f,n=%d,i3=%+d)" % (m, n, i3)


def main():
    N = 8000
    theta, phi, dV = build_grid(N)

    modes = fw.shell_modes(1) + fw.shell_modes(2) + fw.shell_modes(3)
    nm = len(modes)
    E = np.array([n + abs(m) + 0.5 for (m, n, i3) in modes])

    # psi[a] = (2, N) spinor, normalized to unit L2 on the sphere
    psi = np.zeros((nm, 2, N), dtype=complex)
    for a, (m, n, i3) in enumerate(modes):
        u, d = fw.psi(m, n, i3, theta, phi)
        s = np.vstack([u, d])
        nrm = np.sqrt((np.abs(s) ** 2).sum() * dV)
        psi[a] = s / nrm

    S0 = np.array([[1, 0], [0, 1]], dtype=complex)
    S1 = np.array([[0, 1], [1, 0]], dtype=complex)
    S2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    S3 = np.array([[1, 0], [0, -1]], dtype=complex)
    spins = [("s0", S0), ("s1", S1), ("s2", S2), ("s3", S3)]

    # precompute Y_lM on the grid
    Y = {}
    for l in range(3):
        for M in range(-l, l + 1):
            Y[(l, M)] = sph_harm_y(l, M, theta, phi)

    print("# continuum E0+E1 enumeration   N_grid=%d   modes=%d (lambda=1,2,3)" % (N, nm))
    print("# energies in lambda-units: 2E0=2  E0+E1=3  2E1=4")
    print("#")
    for (sname, Sig) in spins:
        Spsi = np.einsum("ij,mjg->mig", Sig, psi)             # Sigma psi_b on grid
        for l in range(3):
            P = np.zeros((nm, nm))
            for M in range(-l, l + 1):
                G = dV * np.einsum("aig,big,g->ab", np.conj(psi), Spsi, Y[(l, M)])
                P += np.abs(G) ** 2
            surv = P > TOL
            if not surv.any():
                print("Sigma=%s  l=%d :  no surviving pairs" % (sname, l))
                continue
            Epair = E[:, None] + E[None, :]
            emin = Epair[surv].min()
            has_2e0 = bool((surv & (np.abs(Epair - 2.0) < 1e-6)).any())
            has_e0e1 = bool((surv & (np.abs(Epair - 3.0) < 1e-6)).any())
            tag = "GROUND=%.0f" % emin
            note = "2E0=%s  E0E1=%s" % ("Y" if has_2e0 else ".", "Y" if has_e0e1 else ".")
            print("Sigma=%s  l=%d :  %s   [%s]" % (sname, l, tag, note))
            # list the pairs at the minimal energy
            aa, bb = np.where(surv & (np.abs(Epair - emin) < 1e-6))
            seen = set()
            for a, b in zip(aa, bb):
                key = (a, b)
                if key in seen:
                    continue
                seen.add(key)
                if P[a, b] < 0.05 * P[surv].max():
                    continue                                   # skip tiny leakage
                print("        E=%.0f  fbar %s  x  f %s   P=%.4f"
                      % (emin, mode_label(*modes[a]), mode_label(*modes[b]), P[a, b]))
        print("#")

    # explicit E0+E1 report for the l=2 scalar (expected clean ground) and effmass
    print("# ---- l=2 scalar (Sigma=s0): E0+E1 contributing pairs and continuum effmass ----")
    Sig = S0
    Spsi = np.einsum("ij,mjg->mig", Sig, psi)
    P = np.zeros((nm, nm))
    for M in range(-2, 3):
        G = dV * np.einsum("aig,big,g->ab", np.conj(psi), Spsi, Y[(2, M)])
        P += np.abs(G) ** 2
    Epair = E[:, None] + E[None, :]
    surv = P > TOL
    tvals = np.arange(0.0, 6.01, 0.5)
    C = np.array([(P[surv] * np.exp(-Epair[surv] * t)).sum() for t in tvals])
    meff = np.log(C[:-1] / C[1:]) / 0.5
    print("#   t      C(t)          m_eff")
    for i in range(len(tvals) - 1):
        print("  %4.1f  %.6e   %6.3f" % (tvals[i], C[i], meff[i]))
    e0e1 = surv & (np.abs(Epair - 3.0) < 1e-6)
    print("#   E0+E1 pairs (l=2 scalar):")
    for a, b in zip(*np.where(e0e1)):
        if P[a, b] < 0.05 * P[e0e1].max():
            continue
        print("        fbar %s  x  f %s   P=%.4f"
              % (mode_label(*modes[a]), mode_label(*modes[b]), P[a, b]))


if __name__ == "__main__":
    main()
