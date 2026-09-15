#!/usr/bin/env python3
# two_meson_gevp_OA_min_free_claude.py  [FREE -- minimal {sigma_00^2, O_A} 2x2 GEVP, NOVAC, no proliferation]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_gevp_OA_min_free_claude.py   (or LREF=2 NVDIR=distill_Nv84)
#
# The clean minimal two-meson resolution: connected (vacuum-subtracted) 2x2 GEVP of
#   sigma_00^2  and  O_A = psibar tilde_tau psi   (mode vertex Phi_00 @ tilde_tau),   rebase t0=4.
# Level 0 -> (2,2) excited single meson (~2 E_1);  level 1 -> two-meson (2 m_sigma).  Proliferation was
# shown redundant (two_meson_prolif_compare_free_claude.py).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = int(os.environ.get("T0", "4"))


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  minimal {sigma_00^2, O_A} 2x2 GEVP (NOVAC, t0=%d)  m_sig~%.3f 2m_sig~%.3f"
          % (tag, dc.L, T0, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    oo = np.outer([o2, oA], [o2, oA])

    C22 = np.zeros(twin); C2A = np.zeros(twin); CAA = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a22 = a2A = aAA = 0.0
        for s in range(ns):
            t = s + dt
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum().real
            a2A += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
        C22[dt], C2A[dt], CAA[dt] = a22 / ns, a2A / ns, aAA / ns
    c = np.zeros((twin, 2, 2))
    for t in range(twin):
        c[t] = np.array([[C22[t], C2A[t] + o2 * oA], [C2A[t] + o2 * oA, CAA[t] + oA ** 2]]) - oo   # connected

    C0 = 0.5 * (c[T0] + c[T0].T)
    wv, Uv = np.linalg.eigh(C0)
    order = np.argsort(wv)[::-1]
    Uk = Uv[:, order] / np.sqrt(wv[order])
    lam = np.full((twin, 2), np.nan)
    for t in range(twin):
        M = Uk.T @ (0.5 * (c[t] + c[t].T)) @ Uk
        lam[t] = np.sort(np.linalg.eigvals(M).real)[::-1]
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])

    print("\n  t  |  level0 (2,2)   level1 (two-meson)")
    for t in range(T0 + 1, twin - 2):
        print("  %2d |   %7.4f        %7.4f" % (t, em[t, 0], em[t, 1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i, (col, mk, lab) in enumerate([("tab:green", "o", "level 0: (2,2)"), ("tab:blue", "s", "level 1: two-meson")]):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=col, marker=mk, ms=4, lw=1.2, label=lab)
    m2 = {1: 0.52, 2: None}.get(dc.L, None)
    for y in [v for v in [m2, 2 * msig] if v is not None]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.axhline(2 * msig, color="tab:blue", ls=":", lw=1.0, alpha=0.6)
    ax.set_ylim(0.3, 1.0)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d minimal {$\sigma^2,O_A$} GEVP: two-meson $\to 2m_\sigma$=%.3f" % (dc.L, 2 * msig))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_OA_min_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
