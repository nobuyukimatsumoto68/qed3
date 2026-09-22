#!/usr/bin/env python3
# split_time_check_claude.py   [FREE FIELD -- splitting-in-time study of the two-meson source]
# Run:  ENS=free python3 split_time_check_claude.py
# Two source sigma's at times (t0, t0+D), two sink sigma's at (t0+T, t0+T+D).  Same 10-diagram contraction
# patterns / coefficients, but the equal-time COLLISION leg tau(s,s) becomes tau(t0,t0+D) -> UNEQUAL time ->
# NO contact (no subtraction needed).  A,B,E only (tadpoles vanish); weights 2*{4,2,2}.  Effmass in T.
# D=0 reproduces the unsplit single-meson; does D>0 lift A+B+E to 2 m_sigma?

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc


def tr(*mats):
    M = mats[0]
    for X in mats[1:]:
        M = M @ X
    return np.trace(M)


def ABE(Phi, tau, a, b, c, d):
    # A = -S_S, B = -T_S, E = C_S^2  with source sigma's at (a,b), sink at (c,d)
    A = -tr(Phi[a], tau[a, b], Phi[b], tau[b, c], Phi[c], tau[c, d], Phi[d], tau[d, a])
    B = -tr(Phi[a], tau[a, c], Phi[c], tau[c, b], Phi[b], tau[b, d], Phi[d], tau[d, a])
    CSa = tr(Phi[a], tau[a, c], Phi[c], tau[c, a])
    CSb = tr(Phi[b], tau[b, d], Phi[d], tau[d, b])
    E = CSa * CSb
    return A, B, E


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE FIELD  splitting-in-time  ncfg=%d" % (tag, len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    styles = {0: ("tab:red", "o"), 1: ("tab:blue", "s"), 2: ("tab:green", "^"), 3: ("tab:purple", "D")}
    for D in [0, 1, 2, 3]:
        Tmax = twin - D - 2
        C = np.zeros(Tmax + 1)
        for T in range(1, Tmax + 1):
            acc = 0.0
            n = 0
            for t0 in range(0, twin - (T + D)):
                a, b, c, d = t0, t0 + D, t0 + T, t0 + T + D
                A, B, E = ABE(Phi, tau, a, b, c, d)
                acc += (2.0 * (4.0 * A + 2.0 * B + 2.0 * E)).real
                n += 1
            C[T] = acc / n
        with np.errstate(all="ignore"):
            em = np.log(C[1:-1] / C[2:])          # em[k] = m_eff at T = k+1
        Ts = np.arange(1, len(em) + 1)
        good = np.isfinite(em) & (C[1:-1] > 0) & (C[2:] > 0)
        col, mk = styles[D]
        ax.plot(Ts[good], em[good], color=col, marker=mk, ms=4, lw=1, label=r"$\Delta=%d$" % D)

        def me(T):
            return em[T - 1] if 0 < T <= len(em) else np.nan
        print("  D=%d:  m_eff(T=6,10,14,18) = %.4f %.4f %.4f %.4f" % (D, me(6), me(10), me(14), me(18)))
    ax.axhline(0.378, color="gray", ls="--", lw=1, alpha=0.6, label=r"$m_\sigma=0.378$")
    ax.axhline(0.756, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2m_\sigma=0.756$")
    ax.set_ylim(0.0, 1.2)
    ax.set_xlabel(r"$T$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("FREE FIELD splitting-in-time: A+B+E, source (0,$\\Delta$) sink (T,T+$\\Delta$)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/split_time_check_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
