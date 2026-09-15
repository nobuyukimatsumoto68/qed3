#!/usr/bin/env python3
# two_meson_two_terms_free_claude.py  [FREE -- PS-PS two-meson: FIRST + SECOND term of Eq (5.5)]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 CONTACT=0.5 python3 two_meson_two_terms_free_claude.py
#
# qed3int_v3-4.pdf Eq (5.1)-(5.6): sigma = eta^dag S xi + xi^dag Stilde eta, PS has S=Stilde=1.
# G4 = <S4> + <Stilde4>.  FIRST term <S4>  = sum_i w_i A_i with legs = D^{-1} = tau.
#      SECOND term <Stilde4> = same 10 patterns with the DH propagator D_ov^{-dag}.
# DH propagator in mode space: tau_DH[a,b] = tau[b,a]^dag  (since D^{-dag}=(D^{-1})^dag; GW: D^{-dag}=1-D^{-1}).
# We compute the TOTAL sum_i w_i A_i for each term (NO extra x2 -- that x2 WAS the second term) and G4 = sum,
# plus diagram A alone.  Contact (CONTACT env) subtracted on equal-time legs of BOTH terms.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

CONTACT = float(os.environ.get("CONTACT", "0.5"))


def main():
    de.CONTACT = CONTACT
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  PS-PS two-meson: 1st (tau) + 2nd (DH=D^-dag) term  CONTACT=%.2f  m_sig~%.3f 2m_sig~%.3f"
          % (tag, dc.L, CONTACT, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    # DH propagator two equivalent forms:
    #  (i)  conj-transpose + time-reverse:  tau_DH[a,b] = tau[b,a]^dagger  (= D^{-dag} = (D^{-1})^dag)
    #  (ii) GW form:  D^{-dag} = 1 - D^{-1}  ->  tau_DH[a,b] = delta_ab I - tau[a,b]
    Nv = tau.shape[-1]
    tau_DH = np.conj(tau).transpose(1, 0, 3, 2)
    tau_DH_gw = -tau.copy()
    for a in range(twin):
        tau_DH_gw[a, a] = tau_DH_gw[a, a] + np.eye(Nv)
    print("# equal-time GW  max|tau[a,a]+tau[a,a]^dag - I| = %.2e" %
          max(np.abs(tau[a, a] + tau[a, a].conj().T - np.eye(Nv)).max() for a in range(twin)))
    print("# DH forms agree  max|tau[b,a]^dag - (delta I - tau)| = %.2e" % np.abs(tau_DH - tau_DH_gw).max())

    def total(leg):
        tot = np.zeros(twin)
        Aonly = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10)
            for s in range(ns):
                acc += de.diags_pair(Phi, leg, s, s + dt).real
            acc /= ns
            tot[dt] = (dc.W10 * acc).sum()
            Aonly[dt] = acc[0]                   # diagram A
        return tot, Aonly

    S4, A1 = total(tau)
    St4, A2 = total(tau_DH)
    G4 = S4 + St4
    Asum = A1 + A2

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))

    print("\n  dt |   S4(1st)     Stilde4(2nd)   G4=sum    |  m_eff: S4     St4     G4   | diagA: A1     A2     A1+A2")
    eS, eSt, eG = eff(S4), eff(St4), eff(G4)
    eA1, eA2, eAs = eff(A1), eff(A2), eff(Asum)
    for dt in range(1, min(24, twin - 1)):
        print("  %2d | %+.3e  %+.3e  %+.3e |  %5.2f  %5.2f  %5.2f  |  %5.2f  %5.2f  %5.2f"
              % (dt, S4[dt], St4[dt], G4[dt], eS[dt], eSt[dt], eG[dt], eA1[dt], eA2[dt], eAs[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(S4, "tab:red", "o", r"1st $\langle S_4\rangle$ ($\tau$)"),
                            (St4, "tab:blue", "s", r"2nd $\langle\tilde S_4\rangle$ (DH)"),
                            (G4, "k", "*", r"$G_4$ = sum")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=4, lw=1.2, label=lab)
    for y in [msig, 2 * msig]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("dt"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d PS-PS: 1st + 2nd term of $G_4$ (contact=%.2f)" % (dc.L, CONTACT))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_two_terms_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
