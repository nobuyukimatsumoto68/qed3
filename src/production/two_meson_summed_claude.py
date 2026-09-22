#!/usr/bin/env python3
# two_meson_summed_claude.py
# Two-meson PS.PS correlator built from the SIGNAL diagrams only: drop F,H,I (noise in the
# DC-subtracted per-diagram plot); keep A,B,C,D,E,G,J with their weights (PS.PS = 2 * W10 . diags).
# Plateau (large-dt constant = vacuum) is fit and subtracted PER JACKKNIFE sample, then the
# subtracted correlator is shown on a LOG plot with jackknife errorbars.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de       # diags_pair (uses CONTACT), LABELS

CONTACT = float(os.environ.get("CONTACT", "0.0"))     # match the informative dcsub plot (raw, contact=0)
DROP = [5, 7, 8]                                       # F, H, I -- pure noise
KEEP = [i for i in range(10) if i not in DROP]        # A,B,C,D,E,G,J
PLAT_LO = int(os.environ.get("PLAT_LO", "24"))        # plateau window [PLAT_LO, twin)


def main():
    de.CONTACT = CONTACT
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  CONTACT=%.2f  keep=%s (drop F,H,I)"
          % (tag, len(dc.KS), CONTACT, [de.LABELS[i] for i in KEEP]))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00

    # per-config summed correlator over the kept diagrams: 2 * sum_{i in keep} W10[i] * D_i(dt)
    wkeep = 2.0 * dc.W10[KEEP]
    C = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        c = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += de.diags_pair(Phi, tau, s, s + dt)
            acc = (acc / ns).real
            c[dt] = float(wkeep @ acc[KEEP])
        C.append(c)
    C = np.array(C)                                   # (ncfg, twin)
    ncfg, twin = C.shape

    # per-jackknife plateau subtraction
    samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])     # (ncfg, twin)
    plat = samp[:, PLAT_LO:].mean(1, keepdims=True)                        # per-sample constant
    Csub = samp - plat
    m = Csub.mean(0)
    err = np.sqrt((ncfg - 1) * np.mean((Csub - m) ** 2, 0))

    # effective mass from the same per-jackknife plateau-subtracted samples
    with np.errstate(all="ignore"):
        em = np.log(Csub[:, :-1] / Csub[:, 1:])       # (ncfg, twin-1)
    emm = em.mean(0)
    eme = np.sqrt((ncfg - 1) * np.mean((em - emm) ** 2, 0))

    print("# plateau window [%d,%d)   <plateau> = %.4e" % (PLAT_LO, twin, plat.mean()))
    print("\n  dt   C_sub          err           |C/err|     m_eff(err)")
    for dt in range(0, twin):
        mestr = "  %7.4f(%.4f)" % (emm[dt], eme[dt]) if dt < twin - 1 else ""
        print("  %2d  %+12.5e  %11.4e   %6.2f%s" % (dt, m[dt], err[dt], abs(m[dt]) / err[dt], mestr))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    pos = m[dts] > 0.0
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.errorbar(dts[pos], m[dts][pos], yerr=err[dts][pos], color="tab:red", marker="o", ms=5,
                lw=1, capsize=2, label="two-meson (A,B,C,D,E,G,J)")
    ax.set_yscale("log")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$C_\mathrm{sub}(dt)$  (plateau-subtracted, per jackknife)")
    ax.set_title("PS two-meson, signal diagrams, plateau-subtracted (%s, %d cfg, contact=%.1f)"
                 % (tag, ncfg, CONTACT))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_summed_c%.1f_%s_claude.png" % (CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)

    # effective-mass plot (valid where C_sub stays positive)
    edt = np.arange(1, twin - 1)
    good = np.isfinite(emm[edt]) & (m[edt] > 0) & (m[edt + 1] > 0)
    at = 0.2                                          # temporal spacing (at0.2 ensembles)
    m_f2_phys = 3.08                                  # 0++ F^2 PHYSICAL mass (glue GEVP, m = acosh/a_t)
    atm_f2 = at * m_f2_phys                            # -> same a_t m lattice units as this effmass plot
    fig2, ax2 = plt.subplots(figsize=(7.5, 5.2))
    ax2.errorbar(edt[good], emm[edt][good], yerr=eme[edt][good], color="tab:red", marker="o", ms=5,
                 lw=1, capsize=2, label="two-meson (A,B,C,D,E,G,J)")
    ax2.axhline(0.35, color="gray", ls="--", lw=1, alpha=0.6, label=r"$a_t m_\sigma\approx0.35$")
    ax2.axhline(0.67, color="gray", ls="-.", lw=1, alpha=0.6, label=r"$2 a_t m_\sigma\approx0.67$")
    ax2.axhline(atm_f2, color="tab:green", ls="-", lw=1.6, alpha=0.85,
                label=r"0++ $F^2$: $a_t m=%.3f$ ($m_\mathrm{phys}=%.2f$)" % (atm_f2, m_f2_phys))
    ax2.set_ylim(0.0, 1.4)
    ax2.set_xlabel(r"$dt$")
    ax2.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    ax2.set_title("PS two-meson effective mass, plateau-subtracted (%s, %d cfg, contact=%.1f)"
                  % (tag, ncfg, CONTACT))
    ax2.legend(fontsize=9)
    fig2.tight_layout()
    out2 = "figs/two_meson_summed_effmass_c%.1f_%s_claude.png" % (CONTACT, tag)
    fig2.savefig(out2, dpi=130)
    plt.close(fig2)
    print("# -> %s" % out2)


if __name__ == "__main__":
    main()
