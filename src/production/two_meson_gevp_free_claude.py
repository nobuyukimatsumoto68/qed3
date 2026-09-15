#!/usr/bin/env python3
# two_meson_gevp_free_claude.py   [FREE FIELD -- {sigma, sigma^2} GEVP with the 1<->2 triangle]
# Run:  ENS=free python3 two_meson_gevp_free_claude.py
# 2x2 GEVP in the basis {sigma (single meson), sigma^2 (two meson)}, contact-subtracted perambulator tilde_tau:
#   C_11 = <sigma(t) sigma(0)>   = -Tr[Phi(0) tau(0,t) Phi(t) tau(t,0)]
#   C_22 = <sigma^2(t) sigma^2(0)> = sum_i W10[i] diag_i(tilde_tau)   (A-J, NM contraction_v3.pdf p.1)
#   C_12 = <sigma(t) sigma^2(0)> = 2 C'  (the triangle; G',I',J' are tadpoles -> 0 under tilde_tau)
#     C' = -Tr[Phi(0) tilde_tau(0,0) Phi(0) tau(0,t) Phi(t) tau(t,0)]     (contraction_v3.pdf p.2)
# Solve C(t) v = lambda C(t0) v; effmass of the two eigenvalues -> should be m_sigma=0.378 and 2m_sigma=0.756.

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
    de.CONTACT = 0.0                                   # raw diags -- subtraction is already in tilde_tau
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE  {sigma, sigma^2} GEVP  t0=%d" % (tag, T0))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv                # tilde_tau: contact removed on equal-time diagonal
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    dS = np.mean([np.trace(Phi[a] @ tt[a, a]).real for a in range(twin)])
    print("# <D_S> after contact subtraction = %.3e (should be ~0)" % dS)

    C11 = np.zeros(twin)
    C12 = np.zeros(twin)
    C22 = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a11 = a12 = a22 = 0.0
        for s in range(ns):
            t = s + dt
            m = Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s]      # meson-propagation block
            a11 += (-np.trace(m)).real                       # C_11 single meson
            a12 += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ m)).real   # 2 C' triangle
            a22 += 2.0 * (dc.W10 * de.diags_pair(Phi, tt, s, t)).sum().real
        C11[dt] = a11 / ns
        C12[dt] = a12 / ns
        C22[dt] = a22 / ns

    # KEEP the vacuum in C_22 (0++): the GEVP resolves it as its own level 0 (E=0) -- no subtraction artifact.
    print("# vacuum KEPT (no subtraction); C_22 large-t plateau = %.4e" % C22[twin - 6:twin - 1].mean())

    # GEVP effmass of the two eigenvalues (t vs t0), then log-ratio in t
    lam = np.full((twin, 2), np.nan)
    C0 = np.array([[C11[T0], C12[T0]], [C12[T0], C22[T0]]])
    for t in range(twin):
        Ct = np.array([[C11[t], C12[t]], [C12[t], C22[t]]])
        try:
            ev = np.sort(eig(Ct, C0)[0].real)[::-1]          # descending
            lam[t] = ev
        except Exception:
            pass
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])                      # (twin-1, 2)

    print("\n  t    lam0      lam1      m_eff0    m_eff1   [m_sig 0.378 ; 2m_sig 0.756]")
    for t in range(T0 + 1, twin - 2):
        print("  %2d   %8.2e %8.2e  %7.4f  %7.4f" % (t, lam[t, 0], lam[t, 1], em[t, 0], em[t, 1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    fig, ax = plt.subplots(figsize=(8, 5.6))
    for i, (lab, col, mk) in enumerate([("GEVP level 0", "tab:red", "o"), ("GEVP level 1", "tab:blue", "s")]):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=col, marker=mk, ms=4, lw=1, label=lab)
    ax.axhline(0.378, color="gray", ls="--", lw=1, alpha=0.6, label=r"$m_\sigma=0.378$")
    ax.axhline(0.756, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2m_\sigma=0.756$")
    ax.set_ylim(0.0, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE {$\sigma,\sigma^2$} GEVP: two levels vs $m_\sigma$, $2m_\sigma$ (t0=%d)" % T0)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_free_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
