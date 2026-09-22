#!/usr/bin/env python3
# t00_wilson_kernel_claude.py  [CHUNK A: the spatial Wilson Hamiltonian kernel W for the fermionic T_00]
#
# W = the FIRST LINE of Eq (IV.1) of the free-limit paper qed3_v2-6.pdf = the spatial Wilson operator on the
# time slice (the Hamiltonian H, NOT the temporal/symplectic term).  Matches dirac_simp.h:269-270 with M5=0:
#   off-diag  W[i,j] = 0.5 kappa_{ij} (-r s0 + gamma(i,j)) exp(i u_{ij}) Omega(i,j)   (j in nn(i))
#   diagonal  W[i,i] = 0.5 r ( sum_{j in nn(i)} kappa_{ij} ) s0                       (the "subtracts 1/2" r-term)
# with gamma(i,j) = cos(alpha_ij) s1 + sin(alpha_ij) s2  (dirac_simp.h:227-230),
#      Omega(i,j) = cos(omega_ij/2) s0 - i sin(omega_ij/2) s3  (dirac_simp.h:232-235, Eq III.1),
#      kappa_{ij} = 2 link_volume_{ij} / ell_{ij} / mean_ell   (dirac_simp.h:360, Eq IV.2), u = 0 (free).
# The M5 term is DROPPED: it is the overlap-projection mass, not the physical energy (NM 2026-09-16).
#
# L1 NOTE: the icosahedron is edge- and vertex-transitive -> kappa is UNIFORM -> it is an overall constant that
# CANCELS in the effmass log-ratio.  So for the free L1 validation we set kappa = 1 (structure only); the exact
# per-link kappa (link_volume/ell) is needed only for L>1 and is left as a TODO.
# See t00_hamiltonian_derivation_claude.md.

import os
import numpy as np
import geom_hopping_claude as gh

S0 = gh.S0
S1 = gh.S1
S2 = gh.S2
S3 = gh.S3


def build_W(geom_dir, L, r=1.0, kappa=None, NS=2, include_diag=True):
    # Returns W (2*nsite, 2*nsite) complex, the spatial Wilson Hamiltonian kernel (u=0, free).
    # include_diag=False omits the on-site +1/2 r kappa 1 diagonal (the pure spin-scalar / sigma-like piece).
    om, alpha, nns, nsite = gh.build(geom_dir, L, NS)
    N = NS * nsite
    W = np.zeros((N, N), complex)
    if kappa is None:                                   # L1 uniform: overall constant, cancels in effmass
        kap = {(i, j): 1.0 for i in range(nsite) for j in nns[i]}
    else:
        kap = kappa
    for i in range(nsite):
        diag = np.zeros((2, 2), complex)
        for j in nns[i]:
            k = kap[(i, j)]
            hop = 0.5 * k * ((-r) * S0 + gh.gamma(alpha, i, j)) @ gh.Omega(om, i, j)
            W[NS * i:NS * i + 2, NS * j:NS * j + 2] += hop
            diag += 0.5 * r * k * S0                     # dirac_simp.h:270 (M5 dropped), per neighbor
        if include_diag:
            W[NS * i:NS * i + 2, NS * i:NS * i + 2] += diag
    return W, nns, nsite


def main():
    GEOM = os.environ.get("GEOM", "../../geometry/data/")
    L = int(os.environ.get("LREF", "1"))
    W, nns, nsite = build_W(GEOM, L)
    N = W.shape[0]
    degs = [len(nns[i]) for i in range(nsite)]
    print("# L=%d  nsite=%d  N=%d  degrees=%s" % (L, nsite, N, sorted(set(degs))))

    # Herm/antiherm split: naive (gamma) part C should be ANTI-hermitian, Wilson (-r s0 + diag) part B hermitian
    # (paper IV.4: C antihermitian -> imaginary spectrum, B hermitian -> real spectrum).
    om, alpha, _, _ = gh.build(GEOM, L)
    C = np.zeros((N, N), complex)                        # naive term: 0.5 kappa gamma Omega
    B = np.zeros((N, N), complex)                        # Wilson term: 0.5 kappa(-r s0)Omega + diagonal
    for i in range(nsite):
        d = np.zeros((2, 2), complex)
        for j in nns[i]:
            C[2 * i:2 * i + 2, 2 * j:2 * j + 2] += 0.5 * (gh.gamma(alpha, i, j)) @ gh.Omega(om, i, j)
            B[2 * i:2 * i + 2, 2 * j:2 * j + 2] += 0.5 * ((-1.0) * S0) @ gh.Omega(om, i, j)
            d += 0.5 * S0
        B[2 * i:2 * i + 2, 2 * i:2 * i + 2] += d
    print("# split check:  W - (C+B) max = %.2e" % np.abs(W - (C + B)).max())
    print("# C anti-herm:  max|C + C^dag| = %.2e   (naive term, expect ~0)" % np.abs(C + C.conj().T).max())
    print("# B hermitian:  max|B - B^dag| = %.2e   (Wilson term, expect ~0)" % np.abs(B - B.conj().T).max())

    # Omega unitarity (already in geom_hopping) + W eigenvalues (real+imag parts)
    ev = np.linalg.eigvals(W)
    print("# eig(W): Re in [% .4f, % .4f]  Im in [% .4f, % .4f]  (B->Re, C->Im)"
          % (ev.real.min(), ev.real.max(), ev.imag.min(), ev.imag.max()))
    print("# W is neither herm nor antiherm (expected): max|W-W^dag|=%.3e" % np.abs(W - W.conj().T).max())


if __name__ == "__main__":
    main()
