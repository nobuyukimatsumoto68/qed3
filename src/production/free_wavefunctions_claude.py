#!/usr/bin/env python3
# free_wavefunctions_claude.py  -- analytic free Dirac eigen-spinors on S^2 (qed3_v2-6.pdf App C.1).
#   xi_{|m|,n}(z) = (1-z)^{(|m|-1/2)/2} (1+z)^{-(|m|+1/2)/2} P^{(|m|-1/2,-|m|-1/2)}_{n+|m|+1/2}(z)   (C.10)
#   psi_{m,n,i3}(th,ph) = e^{i m ph}/sqrt(2pi) * (1/sqrt2) [ xi(im z) ; i3 im i (-1)^n xi(-im z) ]      (C.18)
#   lambda_spatial = n+|m|+1/2 (C.17);  norm (psi,psi) = c^2_{|m|,n}  (C.19-C.20).
# im = sign(m), z = cos(theta).  Verified below against the orthonormality integral.

import numpy as np
from scipy.special import eval_jacobi
from math import gamma, sqrt, pi


def xi(mabs, n, z):
    a = mabs - 0.5
    b = -mabs - 0.5
    deg = int(round(n + mabs + 0.5))
    pref = (1.0 - z) ** (0.5 * (mabs - 0.5)) * (1.0 + z) ** (-0.5 * (mabs + 0.5))
    return pref * eval_jacobi(deg, a, b, z)


def psi(m, n, i3, theta, phi):
    # 2-spinor at (theta,phi); m half-integer, i3=+-1
    mabs = abs(m)
    im = 1.0 if m > 0 else -1.0
    z = np.cos(theta)
    up = xi(mabs, n, im * z)
    dn = i3 * im * 1j * ((-1.0) ** n) * xi(mabs, n, -im * z)
    ph = np.exp(1j * m * phi) / np.sqrt(2.0 * pi) / np.sqrt(2.0)
    return ph * up, ph * dn      # (upper, lower) spinor components


def c2(mabs, n):
    from math import lgamma
    lg = lambda x: lgamma(x + 1.0)
    num = lg(n) + lg(n + 2 * mabs)
    den = lg(n + mabs - 0.5) + lg(n + mabs + 0.5)
    return np.exp(num - den) / (2 * n + 2 * mabs + 1)


def shell_modes(lam):
    # (m, n, i3) with n+|m|+1/2 = lam, both eigenvalue signs i3=+-1
    out = []
    mabs = 0.5
    while mabs <= lam:
        n = lam - mabs - 0.5
        if abs(n - round(n)) < 1e-9 and n >= 0:
            n = int(round(n))
            for msign in (+1, -1):
                for i3 in (+1, -1):
                    out.append((msign * mabs, n, i3))
        mabs += 1.0
    return out


if __name__ == "__main__":
    # verify orthonormality on a Fibonacci sphere grid: (psi,psi') = c^2 delta
    N = 4000
    gi = np.arange(N) + 0.5
    ph = np.mod(2 * pi * gi / ((1 + 5 ** 0.5) / 2), 2 * pi)
    ct = 1 - 2 * gi / N
    th = np.arccos(ct)
    dV = 4 * pi / N                                     # equal-area weight
    modes = shell_modes(1) + shell_modes(2)
    print("# lambda   (m, n, i3)      <psi,psi> (num)     c^2 (analytic)   ratio")
    for (m, n, i3) in modes:
        u, d = psi(m, n, i3, th, ph)
        norm = np.sum((np.abs(u) ** 2 + np.abs(d) ** 2)) * dV
        cc = c2(abs(m), n)
        print("  %.1f     (%+.1f,%d,%+d)     %.6f            %.6f       %.4f"
              % (n + abs(m) + 0.5, m, n, i3, norm.real, cc, norm.real / cc))
    # cross-orthogonality spot check: two different lambda=2 modes
    m1 = modes[4]; m2 = modes[6]
    u1, d1 = psi(*m1, th, ph); u2, d2 = psi(*m2, th, ph)
    ov = np.sum(np.conj(u1) * u2 + np.conj(d1) * d2) * dV
    print("# cross <%s|%s> = %.2e (should be ~0 if m or n differ)" % (m1, m2, abs(ov)))
