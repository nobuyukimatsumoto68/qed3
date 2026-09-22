#!/usr/bin/env python3
# two_meson_diagrams_free_claude.py  [FREE -- diagram-by-diagram two-to-two (sigma^2 -> sigma^2) correlator]
# Run:  ENS=free LREF=2 NVDIR=distill_Nv84 CONTACT=0.5 python3 two_meson_diagrams_free_claude.py
#
# Single-config (free field), translation-averaged, NO jackknife.  The 10 Wick diagrams A-J of
# <sigma_00^2(t) sigma_00^2(0)> (diag_effmass_claude.diags_pair; NM contraction_v3.pdf p.1), improved
# propagator tilde_tau (CONTACT=0.5).  Prints per-diagram correlator + effmass; plots |C_i| (log) and
# the per-diagram effmass with the m_sigma and 2 m_sigma lines.  Weights W10 = {4,2,4,4,2,1,4,1,1,1}.

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
    print("# ENS=%s  FREE L=%d  two-to-two diagram-by-diagram  CONTACT=%.2f  m_sigma~%.3f  2m_sigma~%.3f"
          % (tag, dc.L, CONTACT, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    # per-diagram correlator D[i, dt], translation-averaged over source time.
    # Pass RAW tau: diags_pair applies the equal-time contact itself (de.CONTACT), single subtraction.
    D = np.zeros((10, twin))
    for dt in range(twin):
        ns = twin - dt
        acc = np.zeros(10)
        for s in range(ns):
            acc += de.diags_pair(Phi, tau, s, s + dt).real
        D[:, dt] = acc / ns
    Wsum = 2.0 * np.tensordot(dc.W10, D, axes=(0, 0))    # TOTAL C22 (old x2 norm kept for continuity)

    def effmass(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))

    print("\n  diagram         C(dt=2)      C(dt=6)      C(dt=10)   |  m_eff(4)  m_eff(8)  m_eff(12)")
    for i in range(10):
        e = effmass(D[i])
        print("  %-14s  %+.3e  %+.3e  %+.3e  |  %6.3f   %6.3f   %6.3f"
              % (de.LABELS[i], D[i, 2], D[i, 6], D[i, 10], e[4], e[8], e[12]))
    et = effmass(Wsum)
    print("  %-14s  %+.3e  %+.3e  %+.3e  |  %6.3f   %6.3f   %6.3f"
          % ("TOTAL(2*W.D)", Wsum[2], Wsum[6], Wsum[10], et[4], et[8], et[12]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    # color-blind: distinct marker+color per diagram
    styles = [("tab:red", "o"), ("tab:blue", "s"), ("tab:green", "^"), ("tab:purple", "D"),
              ("tab:orange", "v"), ("tab:brown", "P"), ("tab:pink", "X"), ("tab:olive", "*"),
              ("tab:cyan", "h"), ("k", "+")]

    # panel 1: |C_i(dt)| log
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(10):
        col, mk = styles[i]
        ax.semilogy(dts, np.abs(D[i, 1:]), color=col, marker=mk, ms=3, lw=1, label=de.LABELS[i])
    ax.semilogy(dts, np.abs(Wsum[1:]), color="gray", lw=2.5, alpha=0.7, label="TOTAL")
    ax.set_xlabel("dt")
    ax.set_ylabel(r"$|C_\mathrm{diagram}(dt)|$")
    ax.set_title(r"FREE L=%d two-to-two: per-diagram $|C|$ (contact=%.2f)" % (dc.L, CONTACT))
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out1 = "figs/two_meson_diagrams_corr_free_L%d_claude.png" % dc.L
    fig.savefig(out1, dpi=130)
    plt.close(fig)

    # panel 2: per-diagram effmass (same-sign guard)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(10):
        col, mk = styles[i]
        e = effmass(D[i])
        g = np.isfinite(e) & (D[i, :-1] * D[i, 1:] > 0)
        ax.plot(dts[:-1][g[:len(dts) - 1]] if False else np.arange(1, twin)[:len(e)][g], e[g],
                color=col, marker=mk, ms=3, lw=1, label=de.LABELS[i])
    et = effmass(Wsum)
    gt = np.isfinite(et) & (Wsum[:-1] * Wsum[1:] > 0)
    ax.plot(np.arange(1, twin)[:len(et)][gt], et[gt], color="gray", lw=2.5, alpha=0.8, label="TOTAL")
    for y in [msig, 2 * msig]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("dt")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d two-to-two: per-diagram effmass (contact=%.2f)" % (dc.L, CONTACT))
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    out2 = "figs/two_meson_diagrams_effmass_free_L%d_claude.png" % dc.L
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print("\n# -> %s\n# -> %s" % (out1, out2))


if __name__ == "__main__":
    main()
