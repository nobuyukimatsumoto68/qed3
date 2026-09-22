#!/usr/bin/env python3
# two_meson_vsub_claude.py
# PS two-meson correlator, normal-ordered (contact=1/2) with EXPLICIT vacuum subtraction:
#   C_conn(dt) = <O(t)O(0)> - <O>^2 ,   O = sigma_PS^2 (four-point-consistent one-point),
# vs the plateau subtraction.  O(s) = 2 (D_S,A(s)^2 + D'_S,A(s)) with contact-subtracted equal-time legs
# (S_4 = A^2, S~_4 = B^2, both = same for PS).  Translation-averaged; config jackknife.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de       # S4_pair (with CONTACT), diags_pair

CONTACT = 0.5


def main():
    print("# ENS=%s  ncfg=%d  PS two-meson, explicit <O>^2 vacuum subtraction" % (dc.ENS.split("nu0")[0], len(dc.KS)))
    de.CONTACT = CONTACT
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    FP = []          # four-point 2*G10[tau_NO](dt) per config
    OP = []          # one-point O(s) per config (window mean)
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        Iv = np.eye(tau.shape[-1])
        # four-point (translation-averaged, contact-subtracted)
        fp = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                acc += (dc.W10 * de.diags_pair(Phi, tau, s, s + dt)).sum()
            fp[dt] = 2.0 * (acc / ns).real
        FP.append(fp)
        # one-point O(s) = 2 (D_S,A^2 + D'_S,A) with contact-subtracted equal-time tau
        o = 0.0
        for s in range(twin):
            tss = tau[s, s] - CONTACT * Iv
            DS = np.trace(Phi[s] @ tss)
            DpS = np.trace(Phi[s] @ tss @ Phi[s] @ tss)
            o += (2.0 * (DS ** 2 + DpS)).real
        OP.append(o / twin)
    FP = np.array(FP)          # (ncfg, twin)
    OP = np.array(OP)          # (ncfg,)
    ncfg, twin = FP.shape

    def jk_effmass(Cc):
        n = Cc.shape[0]
        samp = np.array([np.delete(Cc, i, 0).mean(0) for i in range(n)])
        with np.errstate(all="ignore"):
            em = np.log(samp[:, :-1] / samp[:, 1:])
        return samp.mean(0), em.mean(0), np.sqrt((n - 1) * np.mean((em - em.mean(0)) ** 2, 0))

    # (1) plateau subtraction (per config)
    Cplat = FP - FP[:, 24:].mean(1, keepdims=True)
    # (2) explicit <O>^2 subtraction:  C_conn(dt) = <O(t)O(0)> - <O>^2 ; jackknife the whole thing
    n = ncfg
    Cexpl = np.zeros((n, twin))     # jackknife samples of the connected correlator
    for i in range(n):
        mask = np.ones(n, bool)
        mask[i] = False
        fpm = FP[mask].mean(0)
        Obar = OP[mask].mean()
        Cexpl[i] = fpm - 0.5 * Obar ** 2     # <A^2>^2+<B^2>^2 = (1/2)<sigma sigma>^2 (flavor-diagonal vacuum)
    cm_e = Cexpl.mean(0)
    err_e = np.sqrt((n - 1) * np.mean((Cexpl - cm_e) ** 2, 0))
    with np.errstate(all="ignore"):
        em_e = np.log(Cexpl[:, :-1] / Cexpl[:, 1:])
    emm_e = em_e.mean(0)
    eee_e = np.sqrt((n - 1) * np.mean((em_e - emm_e) ** 2, 0))

    cm_p, em_p, ee_p = jk_effmass(Cplat)
    print("  <O> = %.4f (+-%.4f)   <O>^2 = %.4f   (four-point plateau ~ %.4f)"
          % (OP.mean(), OP.std(ddof=1) / np.sqrt(n), OP.mean() ** 2, FP[:, 24:].mean()))
    print("\n  dt   C_expl        |C/err|   m_eff(expl)      m_eff(plateau)")
    for dt in range(1, 18):
        print("  %2d  %+11.4e  %6.2f   %7.4f(%.4f)   %7.4f(%.4f)"
              % (dt, cm_e[dt], abs(cm_e[dt]) / err_e[dt], emm_e[dt], eee_e[dt], em_p[dt], ee_p[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tag = dc.ENS.split("nu0")[0]
    dts = np.arange(1, 16)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(dts, emm_e[dts], yerr=eee_e[dts], color="tab:red", marker="o", ms=5, lw=1, capsize=2,
                label=r"explicit $\langle\sigma\sigma\rangle^2$ subtraction")
    ax.errorbar(dts + 0.12, em_p[dts], yerr=ee_p[dts], color="tab:blue", marker="s", ms=4, lw=1, capsize=2,
                alpha=0.6, label="plateau subtraction")
    ax.axhline(0.67, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2 m_\sigma\approx0.67$")
    ax.axhline(0.35, color="gray", ls="--", lw=1, alpha=0.5, label=r"$m_\sigma\approx0.35$")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    ax.set_ylim(0.0, 1.2)
    ax.set_title("PS two-meson effmass, normal-ordered + <O>^2 vsub (%s, %d cfg)" % (tag, ncfg))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_vsub_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
