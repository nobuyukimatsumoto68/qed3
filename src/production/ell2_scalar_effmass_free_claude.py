#!/usr/bin/env python3
# ell2_scalar_effmass_free_claude.py
#   Lattice free-limit test of the l=2 scalar density  O_2M = sum_x A(x) Y_2M(x) psibar(x) psi(x).
#   Continuum enumeration (e0e1_modes_continuum_claude.py) showed the l=2 scalar (Sigma=identity)
#   has its LOWEST state at E0+E1 (lambda=1 x lambda=2): the 2E0 ground is FORBIDDEN because
#   j=1/2 (x) j=1/2 tops out at l=1 < 2.  Sigma=identity => frame-INDEPENDENT, so no frame transport
#   R(x) is needed -- only the spatial harmonic Y_2M(xhat) at each site.
#
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 ell2_scalar_effmass_free_claude.py
#         ENS=free LREF=2 NVDIR=distill_Nv84 python3 ell2_scalar_effmass_free_claude.py
#
#   Vertex (matched to distill_contract Phi_00, Y00 -> real Y_2M, spin-diagonal):
#     w_2M[j] = dual_areas[j//2] * Y_2M(xhat_{j//2}) ,  Phi_2M(t) = V(t)^dag diag(w_2M) V(t).
#   l=2 channel two-point (M-summed = icosahedral H irrep = rotation-invariant l=2):
#     C_l2(dt) = sum_M  mean_s [ -Tr( Phi_2M(t) tau(t,s) Phi_2M(s) tau(s,t) ) ],  t=s+dt.
#   Reference: sigma_00 (l=0) two-point -> 2E0 = m_sigma.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import math
import numpy as np
import distill_contract_claude as dc


def real_y2m(xhat):
    # 5 real l=2 spherical harmonics on unit vectors xhat (N,3); columns M=-2,-1,0,+1,+2.
    x = xhat[:, 0]
    y = xhat[:, 1]
    z = xhat[:, 2]
    c1 = 0.25 * math.sqrt(5.0 / math.pi)
    c2 = 0.5 * math.sqrt(15.0 / math.pi)
    Y = np.zeros((xhat.shape[0], 5))
    Y[:, 0] = c2 * x * y                       # M=-2
    Y[:, 1] = c2 * y * z                       # M=-1
    Y[:, 2] = c1 * (3.0 * z * z - 1.0)         # M= 0
    Y[:, 3] = c2 * x * z                       # M=+1
    Y[:, 4] = 0.5 * c2 * (x * x - y * y)       # M=+2
    return Y


def twopt(Psnk, Psrc, tau, twin):
    # connected single-meson two-point, source-time averaged:  -Tr[P(t) tau(t,s) P(s) tau(s,t)].
    C = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        acc = 0.0
        for s in range(ns):
            t = s + dt
            acc += (-np.trace(Psnk[t] @ tau[t, s] @ Psrc[s] @ tau[s, t])).real
        C[dt] = acc / ns
    return C


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    e0 = 0.5 * msig
    e0e1 = e0 + 0.26 if dc.L == 1 else np.nan          # lattice E1~0.26 measured at L1 (memory)
    print("# ENS=%s  FREE L=%d  l=2 scalar density (E0+E1)  2E0=%.3f" % (tag, dc.L, msig))

    dual = dc.dual_areas_from_mesh()
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    sites = sites / np.linalg.norm(sites, axis=1, keepdims=True)
    Y2 = real_y2m(sites)                                # (N,5)
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    wl2 = [np.repeat(dual * Y2[:, M], dc.NS) for M in range(5)]

    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Vt = [V[tsrc0 + a].T for a in range(twin)]          # (2Ns, Nv) per timeslice
    Phi00 = [Vt[a].conj().T @ (w00[:, None] * Vt[a]) for a in range(twin)]
    PhiM = [[Vt[a].conj().T @ (wl2[M][:, None] * Vt[a]) for a in range(twin)] for M in range(5)]

    # one-point check: <O_2M> must vanish (Y_2M orthogonal to the uniform loop)
    o1 = [np.mean([(-np.trace(PhiM[M][a] @ (tau[a, a]))).real for a in range(twin)]) for M in range(5)]
    print("# one-point <O_2M> (M=-2..2, should be ~0): %s"
          % np.array2string(np.array(o1), precision=2))

    C00 = twopt(Phi00, Phi00, tau, twin)               # l=0 reference -> 2E0
    Cl2 = np.zeros(twin)
    for M in range(5):
        Cl2 += twopt(PhiM[M], PhiM[M], tau, twin)       # l=2 channel = sum over M

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))
    e2, e0m = eff(Cl2), eff(C00)
    print("\n  dt |  C_l2         m_eff(l=2)  |  m_eff(sigma_00, l=0)")
    for dt in range(1, twin - 1):
        print("  %2d |  %+.4e     %6.3f     |   %6.3f" % (dt, Cl2[dt], e2[dt], e0m[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(Cl2, "tab:red", "o", r"$\ell=2$ scalar  ($E_0+E_1$)"),
                            (C00, "tab:gray", "s", r"$\sigma_{00}$ ($\ell=0$) $=2E_0$")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=4, lw=1.2, label=lab)
    ax.axhline(msig, color="tab:gray", ls="--", lw=1.0, alpha=0.7, label=r"$2E_0=%.3f$" % msig)
    if np.isfinite(e0e1):
        ax.axhline(e0e1, color="tab:red", ls=":", lw=1.0, alpha=0.7, label=r"$E_0+E_1\approx%.3f$" % e0e1)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlabel("dt")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  $\ell=2$ scalar density: single meson at $E_0+E_1$" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/ell2_scalar_effmass_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
