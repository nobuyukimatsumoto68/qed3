#!/usr/bin/env python3
# o_a_twopoint_effmass_free_claude.py  -- effmass of the pure O_A two-point <O_A O_A> (the (2,2) meson).
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 o_a_twopoint_effmass_free_claude.py  (or LREF=2 NVDIR=distill_Nv84)
# O_A = psibar tilde_tau psi (mode vertex Phi_00 @ tilde_tau).  Plots <O_A O_A> (=(2,2)=2E_1) and, for
# contrast, <sigma_00 sigma_00> (=(1,1)=2E_0), with the m_sigma and 2m_sigma reference lines.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  pure O_A two-point effmass  m_sig(2E0)~%.3f" % (tag, dc.L, msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]

    def twopt(P):
        C = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                acc += (-np.trace(P[t] @ tau[t, s] @ P[s] @ tau[s, t])).real
            C[dt] = acc / ns
        return C

    Cs = twopt(Phi)     # <sigma_00 sigma_00>  = (1,1)
    Ca = twopt(PA)      # <O_A O_A>            = (2,2)

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))
    es, ea = eff(Cs), eff(Ca)
    print("\n  dt |  <O_A O_A> corr    m_eff (2,2) | m_eff sigma_00 (1,1)")
    for dt in range(1, twin - 1):
        print("  %2d |  %+.4e      %6.3f       |   %6.3f" % (dt, Ca[dt], ea[dt], es[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(Ca, "tab:green", "o", r"$\langle O_A O_A\rangle$  = (2,2)"),
                            (Cs, "tab:gray", "s", r"$\langle\sigma_{00}\sigma_{00}\rangle$ = (1,1)")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=4, lw=1.2, label=lab)
    ax.axhline(msig, color="tab:gray", ls="--", lw=1.0, alpha=0.7, label=r"$2E_0=m_\sigma=%.3f$" % msig)
    ax.axhline(2 * msig, color="k", ls=":", lw=0.9, alpha=0.5, label=r"$2m_\sigma=%.3f$" % (2 * msig))
    ax.set_ylim(0.2, 1.1)
    ax.set_xlabel("dt"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  pure $O_A$ two-point: $(2,2)$ excited meson" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/o_a_twopoint_effmass_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
