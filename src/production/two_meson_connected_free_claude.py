#!/usr/bin/env python3
# two_meson_connected_free_claude.py   [copy&edit of two_meson_connected_claude.py for the FREE theory]
# Run:  ENS=free python3 two_meson_connected_free_claude.py
# The SAME connected four-point cumulant as yesterday, but on the free field (1 deterministic config, so we
# can go to large t with no noise):
#   C_conn(t) = <O(t)O(0)> - (1/2)<O>^2 - 2 <C_S(t)>^2 ,   sigma = PS - 1/2 (l=0 sum), tilde_tau (CONTACT=0.5).
# Question: at large t, does C_conn plateau at 2 m_sigma, or fall to the single-meson m_sigma (via A)?

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

CONTACT = 0.5


def main():
    de.CONTACT = CONTACT
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE  connected cumulant (l=0, sigma=PS-1/2)  ncfg=%d" % (tag, len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    FP = np.zeros(twin)          # <O(t)O(0)> full four-point (2 sum W10 diags)
    CSt = np.zeros(twin)         # single-sigma two-point C_S(t)
    for dt in range(twin):
        ns = twin - dt
        afp = 0.0
        acs = 0.0
        for s in range(ns):
            t = s + dt
            afp += 2.0 * (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum()
            M = Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s]
            acs += np.trace(M)
        FP[dt] = (afp / ns).real
        CSt[dt] = (acs / ns).real
    # composite one-point <O> = 2(D_S^2 + D'_S) with contact-subtracted equal-time, window mean
    O = 0.0
    for s in range(twin):
        tss = tau[s, s] - CONTACT * Iv
        O += (2.0 * (np.trace(Phi[s] @ tss) ** 2 + np.trace(Phi[s] @ tss @ Phi[s] @ tss))).real
    O /= twin

    Cvac = FP - 0.5 * O ** 2                 # vacuum removed only
    Cconn = FP - 0.5 * O ** 2 - 2.0 * CSt ** 2   # full connected cumulant

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(C[:-1] / C[1:])
    eSingle = eff(CSt)
    eVac = eff(Cvac)
    eConn = eff(Cconn)
    print("  <O> = %.6f   FP(inf~dt28) = %.4e   (1/2)<O>^2 = %.4e" % (O, FP[28], 0.5 * O ** 2))
    print("\n  t   C_S(single)   C_conn (FP-vac-2CS^2)   FP-vac only")
    for t in range(2, twin - 1):
        print("  %2d   %7.4f      %8.4f              %8.4f" % (t, eSingle[t], eConn[t], eVac[t]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, lab, col, mk in [(CSt, r"$C_S$ single", (0.4, 0.4, 0.4), "^"),
                            (Cvac, r"FP $-$ (1/2)$\langle O\rangle^2$ (vac only)", "tab:blue", "s"),
                            (Cconn, r"connected cumulant", "tab:red", "o")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(ts[g], e[g], color=col, marker=mk, ms=4, lw=1, label=lab)
    ax.axhline(0.378, color="gray", ls="--", lw=1, alpha=0.6, label=r"$m_\sigma=0.378$")
    ax.axhline(0.756, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2m_\sigma=0.756$")
    ax.set_ylim(0.0, 1.6)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("FREE connected cumulant (l=0, yesterday's code): large-t behavior")
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_connected_free_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
