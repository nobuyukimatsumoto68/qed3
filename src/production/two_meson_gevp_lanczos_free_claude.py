#!/usr/bin/env python3
# two_meson_gevp_lanczos_free_claude.py   [FREE FIELD -- polyiterated (Lanczos-inflated) {sigma,sigma^2} GEVP]
# Run:  ENS=free python3 two_meson_gevp_lanczos_free_claude.py
# NM asymm_gevp.pdf: inflate the basis with the once-time-shifted operator O_i T -> block-Hankel correlator
#   C_inf(t) = [[C(t), C(t+1), ...], [C(t+1), C(t+2), ...], ...]   (K time-shifts, 2K x 2K),
# each block the 2x2 {sigma, sigma^2} matrix (vacuum KEPT).  Solve C_inf(t) v = lam C_inf(t0) v; the 2K
# eigenvalues -> 2K levels (vacuum, single-ground, single-excited, two-meson, ...).

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
K = 2                                                  # time-shifts -> 2K x 2K inflated basis


def main():
    de.CONTACT = 0.0
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE  polyiterated {sigma,sigma^2} GEVP  K=%d (%dx%d)  t0=%d" % (tag, K, 2 * K, 2 * K, T0))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv                # tilde_tau
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    C11 = np.zeros(twin)
    C12 = np.zeros(twin)
    C22 = np.zeros(twin)
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
    # 2x2 blocks (vacuum kept in C22)
    Cmat = np.zeros((twin, 2, 2))
    Cmat[:, 0, 0] = C11
    Cmat[:, 0, 1] = C12
    Cmat[:, 1, 0] = C12
    Cmat[:, 1, 1] = C22

    def Cinf(t):
        M = np.zeros((2 * K, 2 * K))
        for a in range(K):
            for c in range(K):
                if t + a + c < twin:
                    M[2 * a:2 * a + 2, 2 * c:2 * c + 2] = Cmat[t + a + c]
        return M

    nlev = 2 * K
    lam = np.full((twin, nlev), np.nan)
    C0 = Cinf(T0)
    for t in range(twin - 2 * K):
        try:
            ev = eig(Cinf(t), C0)[0].real
            lam[t] = np.sort(ev)[::-1]                  # descending: largest lam = lowest E
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])

    print("\n  t   " + "  ".join("m%d" % i for i in range(nlev)) + "     [0 ; 0.378 ; ~0.59 ; 0.756]")
    for t in range(T0 + 1, twin - 2 * K - 1):
        print("  %2d  " % t + "  ".join("%6.3f" % em[t, i] if np.isfinite(em[t, i]) else "  nan " for i in range(nlev)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue", "tab:green", "tab:purple", "tab:orange"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % len(cols)], marker="o", ms=3, lw=1, label="level %d" % i)
    for y, lab in [(0.0, "vac"), (0.378, r"$m_\sigma$"), (0.756, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE polyiterated {$\sigma,\sigma^2$} GEVP (K=%d, %dx%d)" % (K, 2 * K, 2 * K))
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_lanczos_free_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
