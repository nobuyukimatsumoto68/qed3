#!/usr/bin/env python3
# two_meson_tildetau_claude.py   [copy&edit of two_meson_cnumber / normord scripts -- NEW IDEA, experimental]
# NM contact-subtraction derivation (contact_subtraction.pdf): the -1/2 in sigma=psibar S psi - 1/2 is a
# per-PROPAGATOR contact, absorbed into a subtracted perambulator
#   tilde_tau(s,s) = V^dag [D_ov^{-1} - 1/2 delta] V = tau(s,s) - 1/2 I   (equal time; off-diagonal raw).
# Build tilde_tau ONCE, then run all diagrams RAW on it (no per-diagram CONTACT).  Interacting Nf2 g0.5.
# Reports per-diagram effmass and the A+B+E summed effmass.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

LABELS = de.LABELS
KEEP = [0, 1, 4]                                    # A, B, E
at = 0.2
atm_f2 = at * 3.08


def main():
    tag = dc.ENS.split("nu0")[0]
    de.CONTACT = 0.0                               # RAW diags -- the subtraction is already IN tilde_tau
    print("# ENS=%s  ncfg=%d  tilde_tau (contact-subtracted perambulator) diagrams" % (tag, len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    C_ABE = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Iv = np.eye(tau.shape[-1])
        tt = tau.copy()
        for a in range(twin):
            tt[a, a] = tau[a, a] - 0.5 * Iv        # build tilde_tau: subtract contact on the equal-time diagonal
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += de.diags_pair(Phi, tt, s, s + dt)
            D[:, dt] = (acc / ns).real
        allD.append(D)
        w = 2.0 * dc.W10[KEEP]
        C_ABE.append(float(w[0]) * D[0] + float(w[1]) * D[1] + float(w[2]) * D[4])
    allD = np.array(allD)
    C_ABE = np.array(C_ABE)
    ncfg, _, twin = allD.shape

    def jk_eff(C):                                 # raw log-ratio effmass, config jk
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])
        with np.errstate(all="ignore"):
            em = np.log(samp[:, :-1] / samp[:, 1:])
        return em.mean(0), np.sqrt((ncfg - 1) * np.mean((em - em.mean(0)) ** 2, 0)), samp.mean(0)

    print("\n  diagram        m(dt6)   m(dt10)  m(dt14)  m(dt18)")
    for i in range(10):
        em, ee, cm = jk_eff(allD[:, i, :])
        def g(dt):
            return em[dt] if np.isfinite(em[dt]) else float("nan")
        print("  %-14s  %6.3f   %6.3f   %6.3f   %6.3f" % (LABELS[i], g(6), g(10), g(14), g(18)))

    emA, eeA, cmA = jk_eff(C_ABE)
    print("\n  A+B+E summed (weights 8,4,4):   [m_sig~0.35 ; 2m_sig~0.7 ; 0++ 0.616]")
    print("  dt   m_eff(err)")
    for dt in range(1, 20):
        print("  %2d   %7.4f(%.4f)" % (dt, emA[dt], eeA[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, 22)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    good = np.isfinite(emA[dts]) & (cmA[dts] > 0) & (cmA[dts + 1] > 0)
    ax.errorbar(dts[good], emA[dts][good], yerr=eeA[dts][good], color="tab:red", marker="o", ms=5, lw=1,
                capsize=2, label=r"A+B+E on $\tilde\tau$")
    styles = {0: ("tab:orange", "v"), 1: ("tab:blue", "s"), 4: ("tab:green", "^")}
    for i in KEEP:
        em, ee, cm = jk_eff(allD[:, i, :])
        g2 = np.isfinite(em[dts])
        col, mk = styles[i]
        ax.errorbar(dts[g2], em[dts][g2], yerr=ee[dts][g2], color=col, marker=mk, ms=3, lw=0.7, alpha=0.6,
                    capsize=1.5, label=LABELS[i])
    ax.axhline(atm_f2, color="tab:green", ls="-", lw=1.4, alpha=0.7, label=r"0++ 0.616")
    ax.axhline(0.7, color="gray", ls="-.", lw=1, alpha=0.5, label=r"$2m_\sigma$ 0.7")
    ax.axhline(0.35, color="gray", ls="--", lw=1, alpha=0.5, label=r"$m_\sigma$ 0.35")
    ax.set_ylim(0.0, 1.4)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"$\tilde\tau$ (contact-subtracted propagator): A+B+E, interacting (%s, %d cfg)" % (tag, ncfg))
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_tildetau_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
