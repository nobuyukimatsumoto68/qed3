#!/usr/bin/env python3
# two_meson_cnumber_ABE_sum_claude.py
# VALIDATION: sum only A(-S_S), B(-T_S), E(C_S^2) with weights 2*W10 = {8,4,4}, c-number diags,
# translation-averaged.  These are connected (decay to 0) -> raw effmass = log C(dt)/C(dt+1), config jk.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import two_meson_cnumber_claude as cn

KEEP = [0, 1, 4]                                  # A, B, E
at = 0.2
atm_f2 = at * 3.08


def main():
    tag = dc.ENS.split("nu0")[0]
    w = 2.0 * dc.W10[KEEP]                        # {8,4,4}
    print("# ENS=%s  ncfg=%d  A+B+E sum, weights %s" % (tag, len(dc.KS), list(w)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    C = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        c = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += cn.diags_pair_cnumber(Phi, tau, s, s + dt)
            acc = (acc / ns).real
            c[dt] = float(w @ acc[KEEP])
        C.append(c)
    C = np.array(C)
    ncfg, twin = C.shape

    samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])
    m = samp.mean(0)
    err = np.sqrt((ncfg - 1) * np.mean((samp - m) ** 2, 0))
    with np.errstate(all="ignore"):
        em = np.log(samp[:, :-1] / samp[:, 1:])
    emm = em.mean(0)
    eme = np.sqrt((ncfg - 1) * np.mean((em - emm) ** 2, 0))

    print("\n  dt   C(dt)          |C/err|   m_eff(err)   [0++ a_t m=0.616 ; 2m_sigma~0.67]")
    for dt in range(0, 20):
        mestr = "  %7.4f(%.4f)" % (emm[dt], eme[dt]) if dt < twin - 1 else ""
        print("  %2d  %+11.4e  %8.2f%s" % (dt, m[dt], abs(m[dt]) / err[dt], mestr))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    edt = np.arange(1, twin - 1)
    good = np.isfinite(emm[edt]) & (m[edt] > 0) & (m[edt + 1] > 0)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.errorbar(edt[good], emm[edt][good], yerr=eme[edt][good], color="tab:red", marker="o", ms=5, lw=1,
                capsize=2, label="A+B+E sum (weights 8,4,4)")
    ax.axhline(atm_f2, color="tab:green", ls="-", lw=1.6, alpha=0.85, label=r"0++ $F^2$: $a_t m=%.3f$" % atm_f2)
    ax.axhline(0.67, color="gray", ls="-.", lw=1, alpha=0.5, label=r"$2 a_t m_\sigma\approx0.67$")
    ax.axhline(0.35, color="gray", ls="--", lw=1, alpha=0.5, label=r"$m_\sigma\approx0.35$")
    ax.set_ylim(0.0, 1.4)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("A+B+E validation effmass, c-number diags (%s, %d cfg)" % (tag, ncfg))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_cnumber_ABE_sum_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
