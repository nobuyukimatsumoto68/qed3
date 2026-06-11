#!/usr/bin/env python3
# compare_offsite_invariants_claude.py
#
# Compare the FRAME-INDEPENDENT part of the continuum vs lattice all-to-all D^-1.
#
# For an OFF-site block P(i0 -> j) the spinor frame is an independent SU(2) on each end
# (continuum: global-sigma_3 spin-connection frame; lattice: link/spin-connection frame),
# so P -> U_{i0} P U_j^dagger with DIFFERENT U on each side.  The only quantities invariant
# under independent left/right SU(2) are the SINGULAR VALUES {s1 >= s2} of the 2x2 block.
# Those are what we compare (equivalently |P|_F = sqrt(s1^2+s2^2) and |det P| = s1 s2).
#
# This reaches exactly the off-diagonal, off-site content that the pure-2D eigenbasis sum
# (check_S2_generic_pair_claude.cu) could not, because the 3D continuum all-to-all is
# damping-protected (e^{-lambda|tau|}) and the lattice all-to-all is an independent reference.
#
# Index map (includes/dirac_ext.h:65):  r = Nx*t + NS*site + spin,  NS = 2,  Nx = NS*n_sites.
# Dm_inv is stored FLAT row-major:  Dm_inv[r*N + c].  We read the 2 source rows (site i0, t=0)
# -- contiguous -- and form the 2x2 block to every (site j, time t).
#
# ORDERING CHECK: the lattice site index must match the continuum's pts_n<L>.dat line order.
# We verify this indirectly: continuum off-site singular values must collapse onto a clean
# function of the geodesic angle gamma (isotropy of S^2); if the lattice ones also track the
# SAME gamma, the orderings are consistent and the comparison is valid.
#
# Usage:
#   python3 compare_offsite_invariants_claude.py CONT.h5 LAT.h5 PTS.dat [i0] [t1,t2,...]
# e.g.
#   python3 compare_offsite_invariants_claude.py cont_prop_L4/Dinv.0.h5 \
#       data_free_vmRe0.000000vmIm0.000000/prop_deter_L4/Dinv.0.h5 \
#       ../../geometry/data/pts_n4.dat 0 1,2,4,8
#
# Read-only; single-thread.

import sys
import numpy as np
import h5py


def meta(path):
    with h5py.File(path, "r") as f:
        return int(f["N"][0]), int(f["n_sites"][0]), int(f["Nt"][0])


def source_rows(path, i0, N):
    # rows r = 2*i0 + b at t=0; flat row-major Dm_inv[r*N + c].
    with h5py.File(path, "r") as f:
        re = f["Dm_inv"]["real"]
        im = f["Dm_inv"]["imag"]
        rows = []
        for b in range(2):
            r = 2 * i0 + b
            sl = slice(r * N, (r + 1) * N)
            rows.append(np.asarray(re[sl]) + 1j * np.asarray(im[sl]))
    return rows


def block(rows, j, t, Nx):
    base = Nx * t + 2 * j
    return np.array([[rows[0][base], rows[0][base + 1]],
                     [rows[1][base], rows[1][base + 1]]], dtype=complex)


def svals(Q):
    return np.linalg.svd(Q, compute_uv=False)   # [s1 >= s2] >= 0


def load_pts(path, ns):
    pts = np.loadtxt(path)
    if pts.shape[0] < ns:
        raise SystemExit("pts file has %d < n_sites=%d lines" % (pts.shape[0], ns))
    pts = pts[:ns, :3]
    pts /= np.linalg.norm(pts, axis=1, keepdims=True)
    return pts


def geodesic(pts, i, j):
    c = float(np.clip(np.dot(pts[i], pts[j]), -1.0, 1.0))
    return np.arccos(c)


def main():
    cont, lat, pts_path = sys.argv[1], sys.argv[2], sys.argv[3]
    i0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    ts = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [1, 2, 4, 8]

    Nc, nsc, Ntc = meta(cont)
    Nl, nsl, Ntl = meta(lat)
    if (Nc, nsc, Ntc) != (Nl, nsl, Ntl):
        raise SystemExit("schema mismatch: cont (%d,%d,%d) vs lat (%d,%d,%d)"
                         % (Nc, nsc, Ntc, Nl, nsl, Ntl))
    N, ns, Nt = Nc, nsc, Ntc
    Nx = 2 * ns
    pts = load_pts(pts_path, ns)

    print("# continuum vs lattice off-site frame-invariants (singular values of 2x2 blocks)")
    print("# N=%d  n_sites=%d  Nt=%d   source site i0=%d (t=0)" % (N, ns, Nt, i0))

    cr = source_rows(cont, i0, N)
    lr = source_rows(lat, i0, N)

    # same-site sanity (j=i0): should be pure sigma_3 (off-diag ~ 0) on both, frame-free.
    for t in ts[:1]:
        Qc = block(cr, i0, t, Nx)
        Ql = block(lr, i0, t, Nx)
        offc = abs(Qc[0, 1]) + abs(Qc[1, 0])
        offl = abs(Ql[0, 1]) + abs(Ql[1, 0])
        diagc = abs(Qc[0, 0]) + abs(Qc[1, 1])
        diagl = abs(Ql[0, 0]) + abs(Ql[1, 1])
        print("# same-site (j=i0) t=%d : cont offdiag/diag=%.2e  lat offdiag/diag=%.2e  (pure sigma_3 -> ~0)"
              % (t, offc / (diagc + 1e-300), offl / (diagl + 1e-300)))

    for t in ts:
        rows = []
        for j in range(ns):
            if j == i0:
                continue
            g = geodesic(pts, i0, j)
            sc = svals(block(cr, j, t, Nx))
            sl = svals(block(lr, j, t, Nx))
            rows.append((g, sc[0], sc[1], sl[0], sl[1]))
        rows.sort(key=lambda r: r[0])

        # ratio lat/cont on the larger singular value, over blocks above a magnitude floor
        s1c = np.array([r[1] for r in rows])
        s1l = np.array([r[3] for r in rows])
        floor = 1e-3 * s1c.max()
        keep = s1c > floor
        ratio = s1l[keep] / s1c[keep]
        rmean = ratio.mean() if keep.any() else float("nan")
        rstd = ratio.std() if keep.any() else float("nan")

        print("#")
        print("## t=%d : off-site blocks sorted by geodesic angle gamma" % t)
        print("#  %8s  %12s %12s  %12s %12s  %10s" %
              ("gamma", "s1_cont", "s2_cont", "s1_lat", "s2_lat", "s1_lat/cont"))
        for (g, a, b, c, d) in rows:
            mark = "" if a > floor else "  (below floor)"
            print("   %8.5f  %12.5e %12.5e  %12.5e %12.5e  %10.4f%s"
                  % (g, a, b, c, d, c / (a + 1e-300), mark))
        print("#  ratio s1_lat/cont over blocks above floor: mean=%.4f  std=%.4f  rel.spread=%.3f"
              % (rmean, rstd, rstd / (abs(rmean) + 1e-300)))
        print("#   (a CONSTANT ratio = structures agree up to the overall 1/(a_s a_t) norm;")
        print("#    smooth gamma-dependence on BOTH sides = site orderings are consistent.)")


if __name__ == "__main__":
    main()
