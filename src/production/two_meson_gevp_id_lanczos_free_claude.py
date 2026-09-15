#!/usr/bin/env python3
# two_meson_gevp_id_lanczos_free_claude.py  [FREE -- {1, sigma, sigma^2} + Lanczos on sigma,sigma^2 -> 5x5]
# Run:  ENS=free python3 two_meson_gevp_id_lanczos_free_claude.py
# Basis {1, sigma, sigma^2, sigma^(1), sigma^2^(1)} (identity NOT time-shifted -> 5x5).  Polyiterated relation
# (NM asymm_gevp.pdf): shifting the source or sink operator by one step -> C at t+1; both shifted -> t+2.
# So the 5x5 is built from the base blocks B(t), B(t+1), B(t+2) with the base
#   B(t) = [[v,0,-v],[0,C11(t),C12(t)],[-v,C12(t),C22(t)]]   (identity carries the vacuum v = <O>^2 plateau).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
from scipy.linalg import eig
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = 3


def main():
    de.CONTACT = 0.0
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE  {1,sigma,sigma^2}+Lanczos(sigma,sigma^2) 5x5 GEVP  t0=%d" % (tag, T0))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    C11 = np.zeros(twin); C12 = np.zeros(twin); C22 = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a11 = a12 = a22 = 0.0
        for s in range(ns):
            t = s + dt
            m = Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s]
            a11 += (-np.trace(m)).real
            a12 += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ m)).real
            a22 += 2.0 * (dc.W10 * de.diags_pair(Phi, tt, s, t)).sum().real
        C11[dt], C12[dt], C22[dt] = a11 / ns, a12 / ns, a22 / ns
    v = C22[twin - 6:twin - 1].mean()
    print("# identity vacuum v = %.4e" % v)

    def B(t):
        return np.array([[v, 0.0, -v], [0.0, C11[t], C12[t]], [-v, C12[t], C22[t]]])

    def M(t):
        m = np.zeros((5, 5))
        B0, B1 = B(t), B(t + 1)
        m[:3, :3] = B0                                 # {1, sigma, sigma^2}
        for i in range(3):
            m[i, 3] = m[3, i] = B1[i, 1]               # coupling to sigma^(1) -> t+1, sigma column
            m[i, 4] = m[4, i] = B1[i, 2]               # coupling to sigma^2^(1) -> t+1
        m[3, 3] = C11[t + 2]                            # both shifted -> t+2
        m[3, 4] = m[4, 3] = C12[t + 2]
        m[4, 4] = C22[t + 2]
        return m

    nlev = 5
    lam = np.full((twin, nlev), np.nan)
    C0 = M(T0)
    for t in range(twin - 2):
        try:
            lam[t] = np.sort(eig(M(t), C0)[0].real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])

    print("\n  t    m0      m1      m2      m3      m4    [0;0.378;0.59;0.756]")
    for t in range(T0 + 1, twin - 3):
        row = "  ".join("%6.3f" % em[t, i] if np.isfinite(em[t, i]) else "  nan" for i in range(nlev))
        print("  %2d   %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue", "tab:green", "tab:purple"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, 0.378, 0.756]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE {$1,\sigma,\sigma^2$}+Lanczos 5x5 GEVP")
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_id_lanczos_free_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
