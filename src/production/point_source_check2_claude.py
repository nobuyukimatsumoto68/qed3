#!/usr/bin/env python3
# point_source_check2_claude.py   [FREE FIELD -- point-source two-meson, WITH vs WITHOUT distillation]
# Run:  ENS=free python3 point_source_check2_claude.py
# Point sources at north pole N (max z) and south S (min z).  Two ways to contract:
#   (a) WITH distillation: mode-space traces of point vertices Phi_x = V^dag P_x V and the perambulator tau.
#   (b) WITHOUT distillation: reconstruct the position-space propagator D(x,ta;y,tb) = V(ta,x) tau V(tb,y)^dag
#       (2x2 spin blocks) at the point sites and contract in position/spin space directly.
# Complete basis (Nv=2Ns=24) => the two MUST agree; a mismatch would mean the distillation trace code is wrong.
# Also prints A, B, E separately (A is the suspected single-meson diagram).  Diagrams summed over N<->S sink.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

NS = dc.NS


def point_phi(Vt, site):
    w = np.zeros(Vt.shape[1])
    w[NS * site] = 1.0
    w[NS * site + 1] = 1.0
    return (Vt.T).conj().T @ (w[:, None] * (Vt.T))


def Dblk(V, tau, ta, x, tb, y):
    # 2x2 spin block of D_ov^{-1}(x,ta ; y,tb) reconstructed from the complete-basis perambulator
    Vx = V[ta][:, NS * x:NS * x + 2].T                  # (2, Nv)
    Vy = V[tb][:, NS * y:NS * y + 2].conj()             # (Nv, 2)
    return Vx @ tau[ta, tb] @ Vy


def tr(*mats):
    M = mats[0]
    for X in mats[1:]:
        M = M @ X
    return np.trace(M)


def main():
    tag = dc.ENS.split("nu0")[0]
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    N = int(np.argmax(sites[:, 2]))
    S = int(np.argmin(sites[:, 2]))
    print("# ENS=%s  FREE FIELD  point sources  N=%d S=%d" % (tag, N, S))
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])

    # per-t: A,B,E via distillation traces and via position-space blocks (summed over N<->S sink)
    Adi = np.zeros(twin); Bdi = np.zeros(twin); Edi = np.zeros(twin)
    Apo = np.zeros(twin); Bpo = np.zeros(twin); Epo = np.zeros(twin)
    CS = np.zeros(twin)
    for t in range(1, twin):
        ad = bd = ed = ap = bp = ep = cs = 0.0
        n = 0
        for t0 in range(0, twin - t):
            a, c = t0, t0 + t
            PN0 = point_phi(V[a], N); PS0 = point_phi(V[a], S)
            PNt = point_phi(V[c], N); PSt = point_phi(V[c], S)
            taa, tac, tcc, tca = tau[a, a], tau[a, c], tau[c, c], tau[c, a]
            for (X, Y) in [(N, S), (S, N)]:               # sink order (N,S) and exchange (S,N)
                PX = point_phi(V[c], X); PY = point_phi(V[c], Y)
                # (a) distillation traces
                ad += (-tr(PN0, taa, PS0, tac, PX, tcc, PY, tca)).real
                bd += (-tr(PN0, tac, PX, tca, PS0, tac, PY, tca)).real
                ed += (tr(PN0, tac, PX, tca) * tr(PS0, tac, PY, tca)).real
                # (b) position-space blocks (N,S source ; X,Y sink)
                ap += (-tr(Dblk(V, tau, a, N, a, S), Dblk(V, tau, a, S, c, X),
                           Dblk(V, tau, c, X, c, Y), Dblk(V, tau, c, Y, a, N))).real
                bp += (-tr(Dblk(V, tau, a, N, c, X), Dblk(V, tau, c, X, a, S),
                           Dblk(V, tau, a, S, c, Y), Dblk(V, tau, c, Y, a, N))).real
                ep += (tr(Dblk(V, tau, a, N, c, X), Dblk(V, tau, c, X, a, N)) *
                       tr(Dblk(V, tau, a, S, c, Y), Dblk(V, tau, c, Y, a, S))).real
            cs += (tr(PN0, tac, PNt, tca)).real
            n += 1
        Adi[t], Bdi[t], Edi[t] = ad / n, bd / n, ed / n
        Apo[t], Bpo[t], Epo[t] = ap / n, bp / n, ep / n
        CS[t] = cs / n

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(C[1:-1] / C[2:])

    di = 2.0 * (4.0 * Adi + 2.0 * Bdi + 2.0 * Edi)
    po = 2.0 * (4.0 * Apo + 2.0 * Bpo + 2.0 * Epo)
    print("\n  DISTILLATION vs POSITION-SPACE agreement (A+B+E total):  max|rel diff| = %.2e"
          % np.nanmax(np.abs((di[1:20] - po[1:20]) / (np.abs(di[1:20]) + 1e-30))))
    print("\n  t   C_S(single)  A(-S_S)   B(-T_S)   E(C_S^2)   A+B+E")
    for t in range(2, min(24, twin - 2)):
        print("  %2d   %7.4f    %7.4f   %7.4f   %7.4f   %7.4f"
              % (t, eff(CS)[t - 1], eff(Adi)[t - 1], eff(Bdi)[t - 1], eff(Edi)[t - 1], eff(di)[t - 1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(1, len(eff(di)) + 1)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for C, lab, col, mk in [(CS, r"$C_S$ single", "gray", "^"), (Edi, "E (two-$\\sigma$)", "tab:green", "D"),
                            (Adi, "A ($-S_S$)", "tab:orange", "v"), (di, "A+B+E total", "tab:red", "o")]:
        e = eff(C)
        g = np.isfinite(e) & (C[1:-1] > 0) & (C[2:] > 0)
        ax.plot(ts[g], e[g], color=col, marker=mk, ms=4, lw=1, label=lab)
    ax.axhline(0.378, color="gray", ls="--", lw=1, alpha=0.6, label=r"$m_\sigma=0.378$")
    ax.axhline(0.756, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2m_\sigma=0.756$")
    ax.set_ylim(0.0, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("FREE FIELD point-source N/S two-meson: A, B, E, total")
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/point_source_check2_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
