#!/usr/bin/env python3
# geom_hopping_claude.py  -- spin connection + covariant point-split vertex on the S2 mesh.
# Loads omega_n{L}.dat / alpha_n{L}.dat (dirac_simp.h format: "i j v"), builds the neighbor list and the
# 2x2 spin-connection matrices, and assembles the point-split scalar vertex
#   O_ps = sum_{<ij>} A_ij  psibar_i  Omega_ij U_ij S  psi_j ,   Omega(i,j)=cos(w/2) s0 - i sin(w/2) s3
# (dirac_simp.h:232-235).  U=1 (free).  Summing over BOTH directed edges makes W^ps Hermitian, since
# Omega(j,i)=Omega(i,j)^dagger (omega_ji=-omega_ij).

import numpy as np

S0 = np.eye(2, dtype=complex)
S1 = np.array([[0, 1], [1, 0]], complex)
S2 = np.array([[0, -1j], [1j, 0]], complex)
S3 = np.array([[1, 0], [0, -1]], complex)


def load_links(path):
    # returns dict{(i,j): value}; for omega also fills (j,i) = -value (dirac_simp.h:35-37)
    d = {}
    with open(path) as f:
        for line in f:
            s = line.split()
            if len(s) < 3:
                continue
            i = int(s[0])
            j = int(s[1])
            v = float(s[2])
            d[(i, j)] = v
    return d


def build(geom_dir, L, NS=2):
    omega = load_links(geom_dir + "omega_n%d.dat" % L)
    alpha = load_links(geom_dir + "alpha_n%d.dat" % L)
    # symmetrize omega: omega[j,i] = -omega[i,j]
    om = dict(omega)
    for (i, j), v in omega.items():
        om.setdefault((j, i), -v)
    # neighbor list from the link keys
    nsites = 1 + max(max(i, j) for (i, j) in om)
    nns = [[] for _ in range(nsites)]
    for (i, j) in om:
        if j not in nns[i]:
            nns[i].append(j)
    return om, alpha, nns, nsites


def Omega(om, i, j):
    w = om[(i, j)]
    return np.cos(0.5 * w) * S0 - 1j * np.sin(0.5 * w) * S3


def gamma(alpha, i, j):
    a = alpha[(i, j)]
    return np.cos(a) * S1 + np.sin(a) * S2


def pointsplit_vertex(geom_dir, L, dual, Y00, NS=2, use_gamma=False, r=1.0):
    # W^ps (2Ns x 2Ns): block (i,j) = A_ij * [Omega  or  0.5(-r s0 + gamma) Omega]  for j in nn(i)
    om, alpha, nns, nsites = build(geom_dir, L, NS)
    N = NS * nsites
    W = np.zeros((N, N), complex)
    for i in range(nsites):
        for j in nns[i]:
            A = Y00 * np.sqrt(dual[i] * dual[j])          # ell=0 link measure (symmetric)
            blk = Omega(om, i, j)
            if use_gamma:
                blk = 0.5 * ((-r) * S0 + gamma(alpha, i, j)) @ blk
            W[NS * i:NS * i + 2, NS * j:NS * j + 2] += A * blk
    return W, nns


if __name__ == "__main__":
    import sys
    GEOM = "../../geometry/data/"
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    om, alpha, nns, nsites = build(GEOM, L)
    print("# L=%d  nsites=%d  #directed links=%d  deg[0]=%d" % (L, nsites, len(om), len(nns[0])))
    # Omega unitarity check
    dev = max(np.abs(Omega(om, i, j) @ Omega(om, i, j).conj().T - S0).max()
              for i in range(nsites) for j in nns[i])
    print("# max|Omega Omega^dag - I| = %.2e (unitary)" % dev)
    # Hermiticity of W^ps
    dual = np.ones(nsites)
    W, _ = pointsplit_vertex(GEOM, L, dual, 1.0)
    print("# W^ps shape=%s  max|W - W^dag| = %.2e (should be ~0)" % (W.shape, np.abs(W - W.conj().T).max()))
