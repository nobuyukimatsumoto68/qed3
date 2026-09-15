#!/usr/bin/env python3
# two_meson_normord_ABE_sum_claude.py
# A+B+E sum effmass, TWO operators side by side:
#   (raw)  sigma^2 with sigma=PS-1/2, level-1 only  -> diags_pair_cnumber (raw tau in connected).
#   (:NO:) :sigma^2: fully normal-ordered            -> diags_pair with tau(s,s)-1/2 I in EVERY leg
#          (subtracts the same-time-slice COLLISION contact -> removes the single-sigma admixture).
# Weights 2*W10[{A,B,E}] = {8,4,4}.  Connected -> raw effmass = log C(dt)/C(dt+1), config jk.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import two_meson_cnumber_claude as cn      # diags_pair_cnumber (raw)
import diag_effmass_claude as de           # diags_pair (tau-1/2 I everywhere = collision-subtracted)

KEEP = [0, 1, 4]
at = 0.2
atm_f2 = at * 3.08


def summed(diag_fn):
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    w = 2.0 * dc.W10[KEEP]
    C = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        c = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += diag_fn(Phi, tau, s, s + dt)
            acc = (acc / ns).real
            c[dt] = float(w @ acc[KEEP])
        C.append(c)
    return np.array(C)


def jk_eff(C):
    n, twin = C.shape
    samp = np.array([np.delete(C, i, 0).mean(0) for i in range(n)])
    with np.errstate(all="ignore"):
        em = np.log(samp[:, :-1] / samp[:, 1:])
    return samp.mean(0), em.mean(0), np.sqrt((n - 1) * np.mean((em - em.mean(0)) ** 2, 0))


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  A+B+E:  raw sigma^2  vs  :sigma^2: (collision contact removed)" % (tag, len(dc.KS)))
    de.CONTACT = 0.5
    Craw = summed(cn.diags_pair_cnumber)
    Cno = summed(de.diags_pair)
    mr, emr, eer = jk_eff(Craw)
    mn, emn, een = jk_eff(Cno)
    print("\n  dt   m_eff(raw sigma^2)   m_eff(:sigma^2:)    [m_sig~0.35 ; 2m_sig~0.67 ; 0++ 0.616]")
    for dt in range(1, 20):
        print("  %2d    %7.4f(%.4f)      %7.4f(%.4f)" % (dt, emr[dt], eer[dt], emn[dt], een[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    twin = Craw.shape[1]
    edt = np.arange(1, twin - 1)
    gr = np.isfinite(emr[edt]) & (mr[edt] > 0) & (mr[edt + 1] > 0)
    gn = np.isfinite(emn[edt]) & (mn[edt] > 0) & (mn[edt + 1] > 0)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.errorbar(edt[gr], emr[edt][gr], yerr=eer[edt][gr], color="tab:red", marker="o", ms=5, lw=1, capsize=2,
                label=r"raw $\sigma^2$ (single-$\sigma$ kept)")
    ax.errorbar(edt[gn] + 0.12, emn[edt][gn], yerr=een[edt][gn], color="tab:blue", marker="s", ms=5, lw=1,
                capsize=2, label=r":$\sigma^2$: collision contact removed")
    ax.axhline(atm_f2, color="tab:green", ls="-", lw=1.6, alpha=0.85, label=r"0++ $F^2$: $a_t m=%.3f$" % atm_f2)
    ax.axhline(0.67, color="gray", ls="-.", lw=1, alpha=0.5, label=r"$2 a_t m_\sigma\approx0.67$")
    ax.axhline(0.35, color="gray", ls="--", lw=1, alpha=0.5, label=r"$m_\sigma\approx0.35$")
    ax.set_ylim(0.0, 1.4)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("A+B+E effmass: raw vs collision-contact-removed (%s, %d cfg)" % (tag, len(dc.KS)))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_normord_ABE_sum_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
