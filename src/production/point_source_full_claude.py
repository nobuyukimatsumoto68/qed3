#!/usr/bin/env python3
# point_source_full_claude.py   [FREE FIELD -- FULL point-source 2-meson, counter term IN THE PROPAGATOR]
# Run:  ENS=free python3 point_source_full_claude.py
# Four bilinears sigma at 1=(N,0) 2=(S,0) 3=(N,t) 4=(S,t)  (N,S = poles).  Full four-point = sum over all
# S_4 Wick contractions of the 2x2 block propagator G(i,j)=D_ov^{-1}(x_i,t_i;x_j,t_j) (reconstructed from
# the complete-basis perambulator).  COUNTER TERM IN THE PROPAGATOR: subtract 1/2 I from each DIAGONAL block
# G(i,i) (the on-site self-contraction, where the GW contact lives); off-diagonal (i!=j) untouched.
# This is the full 10-diagram correlator (all exchanges included) with the on-site term removed at the
# propagator level.  Effmass in t.  Reference: single meson C_S = -Tr[G(1,3)G(3,1)].

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from itertools import permutations
import distill_contract_claude as dc

NS = dc.NS
I2 = np.eye(2, dtype=complex)


def Dblk(V, tau, ta, x, tb, y):
    Vx = V[ta][:, NS * x:NS * x + 2].T
    Vy = V[tb][:, NS * y:NS * y + 2].conj()
    return Vx @ tau[ta, tb] @ Vy


def wick4(G):
    # sum over S_4: each perm -> prod over cycles of (-1) Tr[cycle product]
    total = 0.0 + 0j
    for perm in permutations(range(4)):
        visited = [False] * 4
        val = 1.0 + 0j
        for start in range(4):
            if visited[start]:
                continue
            cyc = []
            i = start
            while not visited[i]:
                visited[i] = True
                cyc.append(i)
                i = perm[i]
            prod = None
            for a in cyc:
                blk = G[a][perm[a]]
                prod = blk if prod is None else prod @ blk
            val *= (-1.0) * np.trace(prod)
        total += val
    return total


def main():
    tag = dc.ENS.split("nu0")[0]
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    N = int(np.argmax(sites[:, 2]))
    S = int(np.argmin(sites[:, 2]))
    print("# ENS=%s  FREE FIELD  FULL point-source 2-meson, counter term in propagator  N=%d S=%d" % (tag, N, S))
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])

    # vacuum <O> = equal-time <sigma(N,t')sigma(S,t')> with counter term (translation-averaged)
    Ovac = 0.0
    for t0 in range(twin):
        g = [[Dblk(V, tau, t0, x, t0, y) for (_, y) in [(t0, N), (t0, S)]] for (_, x) in [(t0, N), (t0, S)]]
        g[0][0] = g[0][0] - 0.5 * I2
        g[1][1] = g[1][1] - 0.5 * I2
        # S_2 Wick: tadpole(N)tadpole(S) + (-1)Tr[G(N,S)G(S,N)]
        Ovac += (np.trace(g[0][0]) * np.trace(g[1][1]) - np.trace(g[0][1] @ g[1][0])).real
    Ovac /= twin

    Cfull = np.zeros(twin)
    Cnosub = np.zeros(twin)
    CS = np.zeros(twin)
    tad = 0.0
    ntad = 0
    for t in range(1, twin):
        acc = accn = cs = 0.0
        n = 0
        for t0 in range(0, twin - t):
            a, c = t0, t0 + t
            # bilinears 1=(N,a) 2=(S,a) 3=(N,c) 4=(S,c)
            pts = [(a, N), (a, S), (c, N), (c, S)]
            G = [[Dblk(V, tau, pts[i][0], pts[i][1], pts[j][0], pts[j][1]) for j in range(4)] for i in range(4)]
            # single-meson reference (before counter term)
            cs += (-np.trace(G[0][2] @ G[2][0])).real
            # counter term IN THE PROPAGATOR: subtract 1/2 I on the diagonal
            Gsub = [[G[i][j].copy() for j in range(4)] for i in range(4)]
            for i in range(4):
                Gsub[i][i] = G[i][i] - 0.5 * I2
                if t == 1 and t0 == 0:
                    tad += np.trace(Gsub[i][i]).real
                    ntad += 1
            acc += wick4(Gsub).real
            accn += wick4(G).real          # without counter term, for comparison
            n += 1
        Cfull[t] = acc / n
        Cnosub[t] = accn / n
        CS[t] = cs / n
    print("# vacuum <O> = %.6f   <O>^2 = %.6e" % (Ovac, Ovac ** 2))
    Cconn = Cfull - Ovac ** 2              # subtract the disconnected vacuum
    print("# on-site tadpole after 1/2 subtraction: Tr[G(i,i)-1/2 I] = %.4e (avg over 4 diagonals)"
          % (tad / max(ntad, 1)))

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(C[1:-1] / C[2:])
    e_full = eff(Cfull)
    e_conn = eff(Cconn)
    e_cs = eff(CS)
    print("\n  t   C_S(single)  FULL(counter-term)   FULL - <O>^2 (vac-sub)")
    for t in range(2, min(24, twin - 2)):
        print("  %2d   %7.4f     %8.4f            %8.4f" % (t, e_cs[t - 1], e_full[t - 1], e_conn[t - 1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(1, len(e_full) + 1)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for C, lab, col, mk in [(CS, r"$C_S$ single", "gray", "^"),
                            (Cfull, "FULL, counter-term (vac-dominated)", "tab:blue", "s"),
                            (Cconn, r"FULL $-\langle O\rangle^2$ (vac-sub)", "tab:red", "o")]:
        e = eff(C)
        g = np.isfinite(e) & (C[1:-1] > 0) & (C[2:] > 0)
        ax.plot(ts[g], e[g], color=col, marker=mk, ms=4, lw=1, label=lab)
    ax.axhline(0.378, color="gray", ls="--", lw=1, alpha=0.6, label=r"$m_\sigma=0.378$")
    ax.axhline(0.756, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2m_\sigma=0.756$")
    ax.set_ylim(0.0, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("FREE FIELD full point-source 2-meson (N/S), counter term in propagator")
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/point_source_full_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
