#!/usr/bin/env python3
# sigma_pair_quantum_numbers_claude.py
#   Continuum classification of the scalar density sigma = psibar psi by (Delta_tot, ell_tot, m_tot).
#   For every ordered mode pair (antifermion abar=(m1,n1,i1), fermion a=(m2,n2,i2)) build the LOCAL
#   density rho(x) = psi_abar(x)^dag Gamma psi_a(x) on S^2 from the analytic free wavefunctions
#   (free_wavefunctions_claude.py, qed3_v2-6.pdf App C.1) and project onto Y_{ell m}:
#       c^{ell m}_{abar,a} = int dOmega  Y_{ell m}^*(x)  rho(x) ,   m_tot = m2 - m1 (azimuthal exact).
#   Report per pair: Delta_tot = lambda1 + lambda2 (lambda-units), m_tot, and the ell's with |c|^2 > TOL.
#   Then aggregate per shell-pair (lambda1,lambda2): which ell appear and the state multiplicity.
#   Gamma = identity (the sigma_00 vertex convention).  Set GAMMA=s3 to test psibar sigma3 psi.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np
from scipy.special import sph_harm_y
import free_wavefunctions_claude as fw

TOL = 1.0e-6
LMAX = 6
GAMMA = os.environ.get("GAMMA", "id")            # "id" (scalar psibar psi) or "s3"


def build_grid(N):
    gi = np.arange(N) + 0.5
    phi = np.mod(2 * np.pi * gi / ((1 + 5 ** 0.5) / 2), 2 * np.pi)
    ct = 1 - 2 * gi / N
    theta = np.arccos(ct)
    dV = 4 * np.pi / N
    return theta, phi, dV


def main():
    N = 12000
    theta, phi, dV = build_grid(N)
    Ylm = {}
    for l in range(LMAX + 1):
        for m in range(-l, l + 1):
            Ylm[(l, m)] = sph_harm_y(l, m, theta, phi)

    s3 = np.array([[1, 0], [0, -1]], complex)
    modes = fw.shell_modes(1) + fw.shell_modes(2) + fw.shell_modes(3)
    # unit-normalize each spinor on the sphere
    spin = {}
    for (m, n, i3) in modes:
        u, d = fw.psi(m, n, i3, theta, phi)
        s = np.vstack([u, d])
        s = s / np.sqrt((np.abs(s) ** 2).sum() * dV)
        spin[(m, n, i3)] = s

    def lam(m, n):
        return int(round(n + abs(m) + 0.5))

    # per-pair decomposition
    print("# sigma = psibar psi   Gamma=%s   density Y_lm decomposition (Delta=lambda1+lambda2)" % GAMMA)
    print("# columns: (l1;m1,n1,i1) x (l2;m2,n2,i2)  Delta  m_tot  ell:|c|^2 ...")
    agg = {}                                       # (lam1,lam2) -> {ell: total weight}
    lines = []
    for (m1, n1, i1) in modes:
        for (m2, n2, i2) in modes:
            L1 = lam(m1, n1)
            L2 = lam(m2, n2)
            if L1 > 2 or L2 > 2:                   # per-pair dump: keep shells 1,2 (Delta<=4); shell 3 aggregated only
                pass
            sa = spin[(m1, n1, i1)]                 # antifermion mode (daggered)
            sb = spin[(m2, n2, i2)]                 # fermion mode
            if GAMMA == "s3":
                sb_g = s3 @ sb
            else:
                sb_g = sb
            rho = np.conj(sa[0]) * sb_g[0] + np.conj(sa[1]) * sb_g[1]   # psi_abar^dag Gamma psi_a
            mtot = int(round(m2 - m1))
            cont = []
            for l in range(abs(mtot), LMAX + 1):
                c = dV * np.sum(np.conj(Ylm[(l, mtot)]) * rho)
                w = abs(c) ** 2
                if w > TOL:
                    cont.append((l, w))
                    key = (L1, L2)
                    agg.setdefault(key, {}).setdefault(l, 0.0)
                    agg[key][l] += w
            if cont and L1 <= 2 and L2 <= 2:
                cs = "  ".join("l%d:%.4f" % (l, w) for (l, w) in cont)
                lines.append("  (%d;%+.1f,%d,%+d) x (%d;%+.1f,%d,%+d)  D=%d  m=%+d  %s"
                             % (L1, m1, n1, i1, L2, m2, n2, i2, L1 + L2, mtot, cs))
    for ln in lines:
        print(ln)

    print("\n# ---- aggregate: which ell appear per shell-pair (lambda1,lambda2), Delta, total weight ----")
    for (L1, L2) in sorted(agg.keys()):
        ells = agg[(L1, L2)]
        ellstr = "  ".join("l%d:%.3f" % (l, ells[l]) for l in sorted(ells))
        print("  (lam1=%d, lam2=%d)  Delta=%d :  %s" % (L1, L2, L1 + L2, ellstr))


if __name__ == "__main__":
    main()
