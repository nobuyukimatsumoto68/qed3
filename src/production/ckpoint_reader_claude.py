#!/usr/bin/env python3
# ckpoint_reader_claude.py -- tested reader for the QED3 U(1) gauge configs (ckpoint_lat.<k>), for qed3-79's
#   Omega-rotated free-kernel build.  Provides: (i) per-timeslice SPATIAL link phases theta; (ii) the edge map
#   (link il <-> (site_a,site_b)) in CKPOINT order = QfeLattice links_n<L>.dat (NOT omega_n<L>.dat, which is a
#   DIFFERENT ordering used for the Dirac hopping); (iii) the S^2 discrete divergence + graph-Laplacian stencil for
#   the Poisson theta-solve; (iv) per-slice total flux + mean plaquette sanity from the primary faces face_n<L>.dat.
#
# FORMAT (verified): ckpoint_lat.<k> = raw native doubles, Nt*(n_links+n_sites); per timeslice s a block of
#   (n_links+n_sites): first n_links = SPATIAL link phases theta (U_i=e^{i theta}), next n_sites = temporal.
#   Spatial link il at timeslice s = flat[s*(n_links+n_sites) + il].  (gauge_ext.h idx_sp/idx_tp; GaugeField::read.)
#   L2: n_links=120, n_sites=42, Nt=128 -> 20736 doubles = 165888 bytes.
# EDGE ORDER (verified): ckpoint spatial link il = line il of geometry/data/links_n<L>.dat = QfeLattice links order.
#   directed_theta(theta_slice, a, b) = +theta_il if links[il]=(a,b), -theta_il if (b,a)  (U(b->a)=U(a->b)^dag).
#
# Run (sanity): python3 ckpoint_reader_claude.py  [ENSDIR] [K] [L]  -- prints shape, |theta| range, mean plaquette
#   (must be > 0 = physical => edge order correct), per-slice total flux histogram.

import os
import sys
import numpy as np

GEOM = os.environ.get("GEOM", "../../geometry/data/")


def n_sites(L):
    return 10 * L * L + 2


def load_links(L):
    # links_n<L>.dat: one line "a b" per spatial link, in QfeLattice (= ckpoint) order.  Returns edges (n_links,2).
    edges = []
    with open(GEOM + "links_n%d.dat" % L) as f:
        for line in f:
            s = line.split()
            if len(s) >= 2:
                edges.append((int(s[0]), int(s[1])))
    return np.array(edges, dtype=int)


def load_faces(L):
    # face_n<L>.dat: one line "a b c" per primary triangle.  Returns (n_faces,3).
    faces = []
    with open(GEOM + "face_n%d.dat" % L) as f:
        for line in f:
            s = line.split()
            if len(s) >= 3:
                faces.append((int(s[0]), int(s[1]), int(s[2])))
    return np.array(faces, dtype=int)


def read_spatial_phases(path, L, Nt=128):
    # returns theta_spatial (Nt, n_links): the U(1) link phases per timeslice.
    edges = load_links(L)
    nl = len(edges)
    ns = n_sites(L)
    flat = np.fromfile(path, dtype="<f8")
    assert flat.size == Nt * (nl + ns), "size %d != %d*(%d+%d)" % (flat.size, Nt, nl, ns)
    block = flat.reshape(Nt, nl + ns)
    return block[:, :nl]                      # spatial part; block[:, nl:] is temporal


class Geom:
    # edge/face geometry + directed-theta / divergence / Laplacian helpers, built once per L.
    def __init__(self, L):
        self.L = L
        self.edges = load_links(L)            # (n_links, 2) directed a->b
        self.faces = load_faces(L)            # (n_faces, 3)
        self.ns = n_sites(L)
        self.nl = len(self.edges)
        # directed lookup (a,b) -> (il, sign)
        self.dmap = {}
        for il, (a, b) in enumerate(self.edges):
            self.dmap[(a, b)] = (il, +1.0)
            self.dmap[(b, a)] = (il, -1.0)
        # graph Laplacian (site x site): L = D - A  (for the Poisson solve Lap(theta) = div(a))
        self.Lap = np.zeros((self.ns, self.ns))
        for (a, b) in self.edges:
            self.Lap[a, a] += 1.0
            self.Lap[b, b] += 1.0
            self.Lap[a, b] -= 1.0
            self.Lap[b, a] -= 1.0

    def directed_theta(self, theta_slice, a, b):
        il, sgn = self.dmap[(a, b)]
        return sgn * theta_slice[il]

    def divergence(self, theta_slice):
        # discrete div at each site: sum_{b~a} theta(a->b).  (RHS of the Poisson solve, up to your sign convention.)
        d = np.zeros(self.ns)
        for il, (a, b) in enumerate(self.edges):
            d[a] += theta_slice[il]
            d[b] -= theta_slice[il]
        return d

    def face_flux(self, theta_slice):
        # oriented flux through each primary triangle = sum of directed theta around (a->b->c->a), wrapped to (-pi,pi].
        f = np.empty(len(self.faces))
        for k, (a, b, c) in enumerate(self.faces):
            s = self.directed_theta(theta_slice, a, b) + self.directed_theta(theta_slice, b, c) + self.directed_theta(theta_slice, c, a)
            f[k] = (s + np.pi) % (2 * np.pi) - np.pi
        return f


def _sanity(ensdir, k, L):
    g = Geom(L)
    path = os.path.join(ensdir, "ckpoint_lat.%d" % k)
    th = read_spatial_phases(path, L)
    Nt = th.shape[0]
    print("# %s : theta_spatial shape %s (Nt x n_links), n_sites=%d n_faces=%d" % (path, th.shape, g.ns, len(g.faces)))
    print("# |theta| max = %.4f  (fraction <pi = %.4f)" % (np.abs(th).max(), (np.abs(th) < np.pi).mean()))
    # plaquette sanity per slice (mean cos(face flux)); physical => > 0.  Random edge order => ~0.
    mc = np.array([np.cos(g.face_flux(th[s])).mean() for s in range(Nt)])
    tot = np.array([g.face_flux(th[s]).sum() / (2 * np.pi) for s in range(Nt)])   # total flux / 2pi per slice
    print("# mean plaquette <cos(flux)> over slices = %.4f +- %.4f   (>0 => links_n%d edge order is CORRECT)"
          % (mc.mean(), mc.std(), L))
    print("# total flux/2pi per slice: min %.3f max %.3f mean %.3f  (near integers => monopole sectors)"
          % (tot.min(), tot.max(), tot.mean()))
    nz = np.argmin(np.abs(tot))
    print("# closest-to-zero-flux slice: s=%d  flux/2pi=%.4f  (use for the trivial-flux sandwich check)" % (nz, tot[nz]))


if __name__ == "__main__":
    ensdir = sys.argv[1] if len(sys.argv) > 1 else "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000"
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    L = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    _sanity(ensdir, k, L)
