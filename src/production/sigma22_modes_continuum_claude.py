#!/usr/bin/env python3
# sigma22_modes_continuum_claude.py
#   Pure-continuum: which (m,n,i3) mode combinations make the (2,2)=2E1 excited SCALAR meson
#   (both fermions in the lambda=2 shell, coupled to l=0).  See e0e1_continuum_impl_plan_claude.md.
#
#   Key fact: the l=0 scalar vertex (Sigma=identity, Y_00 constant) is DIAGONAL in the mode basis,
#     Gamma^00_ab = Y00 int psi_a^dag psi_b = Y00 c^2 delta_ab   (orthonormality),
#   so sigma_00 = psibar psi pairs fermion mode a with the antifermion in the SAME mode a.  Hence:
#     - full sigma_00 correlator  C(t)=sum_a |Gamma_aa|^2 e^{-2E_a t}  -> ground 2E0 (a in lambda=1).
#     - (2,2)=2E1 is the DIAGONAL sum over the lambda=2 shell only: sum_{a in lambda2} psibar_a psi_a
#       = the l=0 (trace) singlet of j=3/2 (x) j=3/2.  No frame transport, no CG mixing needed.
#   Energies in lambda-units E_k=n+|m|+1/2: 2E0=2, E0+E1=3, 2E1=4.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
from scipy.special import sph_harm_y
import free_wavefunctions_claude as fw

TOL = 1.0e-8


def build_grid(N):
    gi = np.arange(N) + 0.5
    phi = np.mod(2 * np.pi * gi / ((1 + 5 ** 0.5) / 2), 2 * np.pi)
    ct = 1 - 2 * gi / N
    theta = np.arccos(ct)
    dV = 4 * np.pi / N
    return theta, phi, dV


def lab(m, n, i3):
    return "(m=%+.1f,n=%d,i3=%+d)" % (m, n, i3)


def main():
    N = 8000
    theta, phi, dV = build_grid(N)
    Y00 = sph_harm_y(0, 0, theta, phi)

    modes = fw.shell_modes(1) + fw.shell_modes(2) + fw.shell_modes(3)
    nm = len(modes)
    E = np.array([n + abs(m) + 0.5 for (m, n, i3) in modes])
    lam2 = np.array([abs(round(n + abs(m) + 0.5) - 2) < 1e-6 for (m, n, i3) in modes])

    psi = np.zeros((nm, 2, N), dtype=complex)
    for a, (m, n, i3) in enumerate(modes):
        u, d = fw.psi(m, n, i3, theta, phi)
        s = np.vstack([u, d])
        psi[a] = s / np.sqrt((np.abs(s) ** 2).sum() * dV)

    # l=0 scalar vertex Gamma_ab = dV sum_grid conj(psi_a).psi_b Y00   (Sigma = identity)
    G = dV * np.einsum("aig,big,g->ab", np.conj(psi), psi, Y00)
    offdiag = np.abs(G - np.diag(np.diag(G))).max()
    print("# l=0 scalar (sigma_00) vertex in mode basis:  max|off-diagonal| = %.2e  (=> diagonal)" % offdiag)
    P = np.abs(G) ** 2                                        # channel strength per pair
    print("# => sigma_00 pairs fermion mode a with antifermion in the SAME mode a; energy 2E_a.\n")

    # (2,2) = 2E1 content: the diagonal lambda=2 modes
    print("# (2,2) = 2E1 mode content (diagonal lambda=2 modes, each: fbar_a x f_a):")
    for a in range(nm):
        if lam2[a] and P[a, a] > TOL:
            print("        2E1  fbar %s  x  f %s   P=%.4f" % (lab(*modes[a]), lab(*modes[a]), P[a, a]))
    n2 = int((lam2 & (np.diag(P) > TOL)).sum())
    print("#   -> %d lambda=2 diagonal modes span the (2,2) l=0 singlet (trace over j=3/2 shell).\n" % n2)

    # continuum effmass: full sigma_00 (all diagonal) vs lambda=2-projected (the (2,2) operator)
    diagP = np.diag(P).copy()
    Efull = E
    Eproj = E[lam2]
    Pproj = diagP[lam2]
    tvals = np.arange(0.0, 6.01, 0.5)
    Cfull = np.array([(diagP * np.exp(-2 * Efull * t)).sum() for t in tvals])
    Cproj = np.array([(Pproj * np.exp(-2 * Eproj * t)).sum() for t in tvals])
    mf = np.log(Cfull[:-1] / Cfull[1:]) / 0.5
    mp = np.log(Cproj[:-1] / Cproj[1:]) / 0.5
    print("#   t   m_eff(sigma_00, all shells)   m_eff((2,2)=lambda2-projected)")
    for i in range(len(tvals) - 1):
        print("  %4.1f        %6.3f                    %6.3f" % (tvals[i], mf[i], mp[i]))
    print("#   (continuum lambda-units: sigma_00 -> 2E0=2 ;  (2,2) projected -> 2E1=4)")


if __name__ == "__main__":
    main()
