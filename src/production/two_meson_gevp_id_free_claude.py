#!/usr/bin/env python3
# two_meson_gevp_id_free_claude.py   [FREE FIELD -- {1, sigma, sigma^2} GEVP with the IDENTITY operator]
# Run:  ENS=free python3 two_meson_gevp_id_free_claude.py
# Add the identity operator to absorb the 0++ vacuum as a genuine basis vector (not a subtraction):
#   C_00 = v (identity two-point, scaled to the vacuum), C_01 = <sigma> = 0,
#   C_02 = -v = a<sigma^2> with a=sqrt(v)  (so C_02^2/C_00 = v = C_22 vacuum plateau).
# 3x3 GEVP; level 0 = vacuum (E=0), levels 1,2 = single-meson + one higher.  tilde_tau throughout.

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
    print("# ENS=%s  FREE  {1, sigma, sigma^2} GEVP  t0=%d" % (tag, T0))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
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

    v = C22[twin - 6:twin - 1].mean()                  # C_22 vacuum plateau <O>^2
    print("# identity: C_00 = v = %.4e, C_02 = -v (matches C_22 vacuum)" % v)

    def Cmat(t):
        return np.array([[v, 0.0, -v],
                         [0.0, C11[t], C12[t]],
                         [-v, C12[t], C22[t]]])

    nlev = 3
    lam = np.full((twin, nlev), np.nan)
    C0 = Cmat(T0)
    for t in range(twin):
        try:
            ev = eig(Cmat(t), C0)[0].real
            lam[t] = np.sort(ev)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])

    print("\n  t    m0       m1       m2      [0 ; 0.378 ; 0.59/0.756]")
    for t in range(T0 + 1, twin - 2):
        row = "  ".join("%7.4f" % em[t, i] if np.isfinite(em[t, i]) else "   nan " for i in range(nlev))
        print("  %2d   %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, 0.378, 0.756]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE {$1,\sigma,\sigma^2$} GEVP (identity absorbs the vacuum)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_id_free_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
