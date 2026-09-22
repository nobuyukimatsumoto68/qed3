#!/usr/bin/env python3
# ground_coupling_check_claude.py
#   Does O_A (= diagram-A kernel, psibar tilde_tau psi) couple to the GROUND single meson (2E0)?
#   Measure the cross <sigma_00(t) O_A(0)> and the diagonal <sigma_00(t) sigma_00(0)>, single free config.
#   The large-t ratio R(t) = C_sigmaA / C_sigmasigma -> <0|O_A|ground> / <ground|sigma_00|0> = O_A's
#   ground overlap RELATIVE to sigma_00.  R -> 0 means O_A is orthogonal to the ground (no coupling);
#   R != 0 means O_A couples to the single ground meson.
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 ground_coupling_check_claude.py
#        ENS=free LREF=2 NVDIR=distill_Nv84 python3 ground_coupling_check_claude.py

import os
os.environ.setdefault("ENS", "free")
os.environ.setdefault("LREF", "1")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

LTAG = os.environ.get("LTAG", "L%d" % dc.L)


def build():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]

    Css = np.full(twin, np.nan)
    CsA = np.full(twin, np.nan)
    C12 = np.full(twin, np.nan)                 # triangle <sigma_00(t) sigma^2(0)>: source sigma^2 self-contracts (tt), sink sigma_00
    for dt in range(twin):
        ns = twin - dt
        if ns <= 0:
            continue
        ass = 0.0
        asa = 0.0
        a12 = 0.0
        for s in range(ns):
            t = s + dt
            ass += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            asa += (-np.trace(Phi[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
            a12 += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s])).real
        Css[dt] = ass / ns
        CsA[dt] = asa / ns
        C12[dt] = a12 / ns
    return Css, CsA, C12, twin


def main():
    Css, CsA, C12, twin = build()
    R = CsA / Css
    R12 = C12 / Css
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    with np.errstate(all="ignore"):
        em12 = np.log(np.abs(C12[:-1] / C12[1:]))
    print("# %s: TRIANGLE <sigma_00(t) sigma^2(0)> = C12  (the sigma^2 -> single-meson overlap)  m_sig=%.3f" % (LTAG, msig))
    print("# also C_sA=<sigma_00 O_A>.  C12 -> 0 : sigma^2 does NOT couple to a single meson via this triangle.")
    print("#  dt |   C_ss          C12(triangle)   R12=C12/Css   effmass(C12)   [C_sA/Css]")
    for dt in range(1, twin - 1):
        if np.isfinite(Css[dt]) and abs(Css[dt]) > 0:
            e = em12[dt] if dt < twin - 1 and np.isfinite(em12[dt]) else np.nan
            print("#  %2d | %+.5e   %+.5e   %+.6f   %7.4f   [%+.6f]" % (dt, Css[dt], C12[dt], R12[dt], e, R[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(twin)
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    ax.axhline(0.0, color="gray", lw=1, alpha=0.7)
    g = np.isfinite(R12)
    ax.plot(ts[g], R12[g], color="tab:red", marker="o", ms=5, lw=1.2, label=r"$R_{12}=\langle\sigma_{00}\sigma^2\rangle/\langle\sigma_{00}\sigma_{00}\rangle$ (triangle)")
    g2 = np.isfinite(R)
    ax.plot(ts[g2], R[g2], color="tab:blue", marker="s", ms=4, lw=1.0, label=r"$\langle\sigma_{00}O_A\rangle/\langle\sigma_{00}\sigma_{00}\rangle$")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$R(t)$  ($\sigma^2\to$single-meson coupling)")
    ax.set_title(r"FREE %s: does $O_A$ couple to the ground single meson?" % LTAG)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/ground_coupling_%s_claude.png" % LTAG
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
