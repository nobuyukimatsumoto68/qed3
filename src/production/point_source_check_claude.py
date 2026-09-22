#!/usr/bin/env python3
# point_source_check_claude.py   [FREE FIELD -- point-source two-to-two meson correlator]
# Run:  ENS=free python3 point_source_check_claude.py
# Two POINT sources at the north pole N (max z) and south pole S (min z), NOT the l=0 sum.  Source sigma's
# sigma(N,0) sigma(S,0), sink sigma(N,t) sigma(S,t).  N != S => the equal-time collision tau(N,S) has NO
# contact (delta(N,S)=0), so NO subtraction needed.  Same A,B,E contraction patterns / coefficients, SUMMED
# over the x<->y (N<->S sink) exchange.  Complete basis (Nv=2Ns=24) => exact.  Compare vs single-meson C_S(N->N).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

NS = dc.NS


def point_phi(Vt, site):
    # point vertex at a single site (sum over the 2 spin components), same construction as the l=0 Phi
    w = np.zeros(Vt.shape[1])
    w[NS * site] = 1.0
    w[NS * site + 1] = 1.0
    return (Vt.T).conj().T @ (w[:, None] * (Vt.T))


def tr(*mats):
    M = mats[0]
    for X in mats[1:]:
        M = M @ X
    return np.trace(M)


def ABE(Pa1, Pa2, Pc1, Pc2, taa, tac, tcc, tca):
    A = -tr(Pa1, taa, Pa2, tac, Pc1, tcc, Pc2, tca)
    B = -tr(Pa1, tac, Pc1, tca, Pa2, tac, Pc2, tca)
    E = tr(Pa1, tac, Pc1, tca) * tr(Pa2, tac, Pc2, tca)
    return A, B, E


def main():
    tag = dc.ENS.split("nu0")[0]
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    N = int(np.argmax(sites[:, 2]))
    S = int(np.argmin(sites[:, 2]))
    print("# ENS=%s  FREE FIELD  point sources  N=%d (z=%.3f)  S=%d (z=%.3f)"
          % (tag, N, sites[N, 2], S, sites[S, 2]))
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])

    Ctwo = np.zeros(twin)          # symmetrized 2-meson A+B+E
    Cone = np.zeros(twin)          # single meson C_S(N->N)
    for t in range(1, twin):
        acc = 0.0
        acc1 = 0.0
        nsrc = 0
        for t0 in range(0, twin - t):
            a, c = t0, t0 + t
            PN_s = point_phi(V[a], N)
            PS_s = point_phi(V[a], S)
            PN_t = point_phi(V[c], N)
            PS_t = point_phi(V[c], S)
            taa, tac, tcc, tca = tau[a, a], tau[a, c], tau[c, c], tau[c, a]
            # direct sink (N,S) + exchange sink (S,N)
            for Pc1, Pc2 in [(PN_t, PS_t), (PS_t, PN_t)]:
                A, B, E = ABE(PN_s, PS_s, Pc1, Pc2, taa, tac, tcc, tca)
                acc += (2.0 * (4.0 * A + 2.0 * B + 2.0 * E)).real
            acc1 += (-tr(PN_s, tac, PN_t, tca)).real     # single meson N->N (= -C_S loop)
            nsrc += 1
        Ctwo[t] = acc / nsrc
        Cone[t] = acc1 / nsrc

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(C[1:-1] / C[2:])
    e2 = eff(Ctwo)
    e1 = eff(Cone)
    print("\n  t    C_S(single)   A+B+E(two, N/S point)")
    for t in range(1, min(24, twin - 2)):
        print("  %2d    %7.4f       %7.4f" % (t, e1[t - 1] if t - 1 < len(e1) else np.nan,
                                              e2[t - 1] if t - 1 < len(e2) else np.nan))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(1, len(e2) + 1)
    g2 = np.isfinite(e2) & (Ctwo[1:-1] > 0) & (Ctwo[2:] > 0)
    g1 = np.isfinite(e1) & (Cone[1:-1] > 0) & (Cone[2:] > 0)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(ts[g1], e1[g1], color="gray", marker="^", ms=3, lw=0.8, label=r"$C_S$ single ($m_\sigma$)")
    ax.plot(2 * ts[g1], 2 * e1[g1] * 0 + 2 * e1[g1], alpha=0)  # placeholder
    ax.plot(ts[g2], e2[g2], color="tab:red", marker="o", ms=4, lw=1, label="two-meson N/S point (A+B+E)")
    ax.axhline(0.378, color="gray", ls="--", lw=1, alpha=0.6, label=r"$m_\sigma=0.378$")
    ax.axhline(0.756, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2m_\sigma=0.756$")
    ax.set_ylim(0.0, 1.2)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("FREE FIELD point-source two-meson (N=north, S=south pole)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/point_source_check_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
