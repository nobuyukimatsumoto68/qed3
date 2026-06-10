#!/usr/bin/env python3
# compare_gsgt_samesite_claude.py
#
# Same-site G_s/G_t(t) computed DIRECTLY from a saved all-to-all propagator
# (Dm_inv) over a t range (default 0 < t < 12), as an INDEPENDENT reference for
# the other agent's jj contraction code.
#
# Contraction (vector current, qed3int_v2-11 Sec.4):
#   f^{ab}(t) = -tr[ sigma^a G(x1,x2) sigma^b G(x2,x1) ],  x1=(i0,0), x2=(i0,dt).
#   G_t = f^33 ,  G_s = f^11 + f^22 .  CFT (Delta=2): G_s/G_t = -(D-1) = -2, t-independent.
#
# Both directions G(x1,x2) AND G(x2,x1) are read straight out of the DENSE all-to-all
# matrix -- NO time-translation / antiperiodic-image assumption (that is precisely the
# step a contraction code can get wrong).  The same-site continuum block is pure
# c3(tau) sigma_3 with c3 ODD in tau, so analytically f=diag(-2,-2,+2) c3^2 and
# G_s/G_t = -2 exactly at every t, regardless of c3.
#
# Index map (includes/dirac_ext.h:65):  r = Nx*t + 2*site + spin ,  Nx = 2*n_sites .
# Dm_inv stored FLAT row-major:  Dm_inv[r*N + c].  (row<->col swap only sends f^{ab}->f^{ba},
# so the diagonal G_t, G_s are convention-independent.)
#
# Usage:   python3 compare_gsgt_samesite_claude.py DINV.h5 [i0] [a_t] [t_max]
#   e.g.   python3 compare_gsgt_samesite_claude.py cont_prop_L1/Dinv.0.h5 0 0.2 12.0
#
# Read-only; single-thread.

import sys
import numpy as np
import h5py

SIG = [np.array([[0, 1], [1, 0]], dtype=complex),       # sigma_1
       np.array([[0, -1j], [1j, 0]], dtype=complex),     # sigma_2
       np.array([[1, 0], [0, -1]], dtype=complex)]       # sigma_3


def main():
    path = sys.argv[1]
    i0   = int(sys.argv[2])   if len(sys.argv) > 2 else 0
    at   = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2
    tmax = float(sys.argv[4]) if len(sys.argv) > 4 else 12.0

    with h5py.File(path, "r") as f:
        N  = int(f["N"][0])
        ns = int(f["n_sites"][0])
        Nt = int(f["Nt"][0])
        Nx = 2 * ns
        re = f["Dm_inv"]["real"]
        im = f["Dm_inv"]["imag"]

        # rows for sink (i0, t=0): full length N, contiguous.
        R0 = np.empty((2, N), dtype=complex)
        for b in range(2):
            r = 2 * i0 + b
            R0[b] = np.asarray(re[r * N:(r + 1) * N]) + 1j * np.asarray(im[r * N:(r + 1) * N])

        dt_max = int(round(tmax / at))
        print("# same-site G_s/G_t  from %s" % path)
        print("# N=%d  n_sites=%d  Nt=%d   i0=%d   a_t=%.3f   (CFT: G_s/G_t = -2)" %
              (N, ns, Nt, i0, at))
        print("#   t        G_t(f33)        G_s(f11+f22)     G_s/G_t      dev")
        for dt in range(1, dt_max + 1):
            t = dt * at
            base2 = Nx * dt + 2 * i0

            # A = G(sink=(i0,0), source=(i0,dt))   [2x2]
            A = R0[:, base2:base2 + 2]

            # B = G(sink=(i0,dt), source=(i0,0))   [2x2]
            B = np.empty((2, 2), dtype=complex)
            for b in range(2):
                r = Nx * dt + 2 * i0 + b
                lo = r * N + 2 * i0
                B[b] = np.asarray(re[lo:lo + 2]) + 1j * np.asarray(im[lo:lo + 2])

            f33 = -np.trace(SIG[2] @ A @ SIG[2] @ B)
            f11 = -np.trace(SIG[0] @ A @ SIG[0] @ B)
            f22 = -np.trace(SIG[1] @ A @ SIG[1] @ B)

            Gt = f33.real
            Gs = (f11 + f22).real
            ratio = Gs / Gt if Gt != 0.0 else float("nan")
            print("  %5.2f   %+.7e   %+.7e   %+.6f   %+.1e" %
                  (t, Gt, Gs, ratio, ratio + 2.0))


if __name__ == "__main__":
    main()
