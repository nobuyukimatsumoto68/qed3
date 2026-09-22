#!/usr/bin/env python3
# state1111_freelimit_claude.py
#   Chunk 1: free-limit classification of the {1,1,1,1} level and a rotational-covariance
#   check on reading (B) = "both legs (n=1,|m|=1/2)".
#
#   Free Dirac shell on S^2 (qed3_v2-6.pdf App C, Eq C.17): lambda = n + |m| + 1/2,
#   single-particle j = n + |m| = lambda - 1/2.  A shell lambda has degeneracy 4*lambda.
#   If the shell is ONE j=(lambda-1/2) multiplet (times iota_3=+-1), then 2*(2j+1) = 4*lambda,
#   which it is -- so (n,|m|) are NOT independent: within a shell, n = lambda - 1/2 - |m|.
#
#   Question tested here: is the ell=0 scalar restricted to the (n=1,|m|=1/2) members of the
#   lambda=2 shell a PURE ell=0 state, or does dropping the (n=0,|m|=3/2) members (the |m|=3/2
#   part of the SAME j=3/2 multiplet) leak into ell=2?  Method: build the shell-trace scalar
#   density rho(x) = sum_a psi_a^dag psi_a over the chosen mode set and project onto Y_{ell,0}.

import numpy as np
from math import pi, sqrt
try:
    from scipy.special import sph_harm_y
    def _Ylm(ell, m, th, ph):
        return sph_harm_y(ell, m, th, ph)
except ImportError:
    from scipy.special import sph_harm
    def _Ylm(ell, m, th, ph):
        return sph_harm(m, ell, ph, th)
import free_wavefunctions_claude as fw


def enumerate_shell(lam):
    modes = fw.shell_modes(lam)
    print("# lambda=%d  deg=%d   modes (m, n, i3),  j=n+|m|:" % (lam, len(modes)))
    for (m, n, i3) in modes:
        j = n + abs(m)
        print("    (m=%+.1f, n=%d, i3=%+d)   j=%.1f" % (m, n, i3, j))
    return modes


# Fibonacci sphere grid (equal-area)
N = 20000
gi = np.arange(N) + 0.5
phi = np.mod(2 * pi * gi / ((1 + 5 ** 0.5) / 2), 2 * pi)
ct = 1 - 2 * gi / N
theta = np.arccos(ct)
dV = 4 * pi / N


def shell_trace_density(modeset):
    # rho(x) = sum_{a in set} psi_a^dag(x) psi_a(x)  (Gamma = identity, diagonal a=a)
    rho = np.zeros(N, dtype=complex)
    for (m, n, i3) in modeset:
        u, d = fw.psi(m, n, i3, theta, phi)
        cc = fw.c2(abs(m), n)
        rho += (np.abs(u) ** 2 + np.abs(d) ** 2) / cc      # unit-normalized modes
    return rho


def ell_content(rho, ellmax=4):
    # overlap with Y_{ell,0} (m_tot=0 sector); scipy sph_harm(m, ell, phi, theta)
    out = []
    for ell in range(0, ellmax + 1):
        Y = _Ylm(ell, 0, theta, phi)
        c = np.sum(np.conj(Y) * rho) * dV
        out.append((ell, abs(c)))
    return out


if __name__ == "__main__":
    print("===== shell enumeration =====")
    for lam in (1, 2, 3):
        enumerate_shell(lam)
        print()

    print("===== lambda=2 ell-content of the shell-trace scalar =====")
    shell2 = fw.shell_modes(2)
    full = shell2
    n1_only = [(m, n, i3) for (m, n, i3) in shell2 if n == 1]      # (n=1,|m|=1/2)
    n0_only = [(m, n, i3) for (m, n, i3) in shell2 if n == 0]      # (n=0,|m|=3/2)

    for label, ms in (("FULL lambda=2 shell", full),
                      ("(n=1,|m|=1/2) only  ", n1_only),
                      ("(n=0,|m|=3/2) only  ", n0_only)):
        rho = shell_trace_density(ms)
        cont = ell_content(rho)
        s = "  ".join("ell%d=%.4f" % (e, v) for (e, v) in cont)
        print("%s :  %s" % (label, s))
