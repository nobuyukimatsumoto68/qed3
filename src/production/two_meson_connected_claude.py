#!/usr/bin/env python3
# two_meson_connected_claude.py
# FULLY-CONNECTED four-point cumulant of sigma = PS - 1/2 (Chester-Pufu two-scalar 0++):
#   <sigma sigma sigma sigma>_c = <O(t)O(0)> - <O>^2 - 2 <C_S(t)>^2 ,   O = sigma^2 (= two-meson op).
# The three pairings of <sigma(x)sigma(y)sigma(z)sigma(w)> (x,y at t; z,w at 0):
#   (xy)(zw) = <O(t)><O(0)> = <O>^2         (vacuum / one-point squared)
#   (xz)(yw)+(xw)(yz) = 2 <C_S(t)>^2         (two FREE single-sigma propagators)
# All averages are ENSEMBLE means; the whole combination is config-jackknifed.  sigma = PS - 1/2
# => equal-time tau(s,s) contact-subtracted by 1/2 (CONTACT=0.5); C_S uses only off-diagonal legs.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de       # diags_pair (CONTACT), LABELS

CONTACT = 0.5                          # sigma = PS - 1/2
PLAT_LO = 24
at = 0.2
atm_f2 = at * 3.08


def main():
    de.CONTACT = CONTACT
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  connected four-point cumulant (sigma = PS - 1/2)" % (tag, len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    FP = []           # <O(t)O(0)> full four-point (all 10 diags), per config
    CS = []           # C_S(t) single-sigma two-point, per config
    OP = []           # O one-point (window mean), per config
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Iv = np.eye(tau.shape[-1])
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        fp = np.zeros(twin)
        cs = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            afp = 0.0
            acs = 0.0
            for s in range(ns):
                t = s + dt
                afp += 2.0 * (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum()
                M = Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s]
                acs += np.trace(M)
            fp[dt] = (afp / ns).real
            cs[dt] = (acs / ns).real
        o = 0.0
        for s in range(twin):
            tss = tau[s, s] - CONTACT * Iv
            o += (2.0 * (np.trace(Phi[s] @ tss) ** 2 + np.trace(Phi[s] @ tss @ Phi[s] @ tss))).real
        FP.append(fp)
        CS.append(cs)
        OP.append(o / twin)
    FP = np.array(FP)
    CS = np.array(CS)
    OP = np.array(OP)
    ncfg, twin = FP.shape

    # config jackknife.  Vacuum coeff = 1/2 (S_4+S~_4 flavor factor; FP(inf) = (1/2)<O>^2).
    #   V0 = FP - (1/2)<O>^2                  (vacuum removed only)
    #   V1 = FP - (1/2)<O>^2 - 2 <C_S>^2      (full connected cumulant: also free two-sigma)
    def jk_combo(free_coeff):
        s = np.zeros((ncfg, twin))
        for i in range(ncfg):
            m = np.ones(ncfg, bool)
            m[i] = False
            s[i] = FP[m].mean(0) - 0.5 * OP[m].mean() ** 2 - free_coeff * CS[m].mean(0) ** 2
        c = s.mean(0)
        e = np.sqrt((ncfg - 1) * np.mean((s - c) ** 2, 0))
        with np.errstate(all="ignore"):
            em = np.log(s[:, :-1] / s[:, 1:])
        return c, e, em.mean(0), np.sqrt((ncfg - 1) * np.mean((em - em.mean(0)) ** 2, 0))

    Cc0, e0, em0, ee0 = jk_combo(0.0)          # vacuum-only
    Cc, err, emm, eme = jk_combo(2.0)          # full cumulant
    print("  <O> = %.5f   (1/2)<O>^2 = %.4e   FP(inf~dt30) = %.4e" % (OP.mean(), 0.5 * OP.mean() ** 2, FP.mean(0)[30]))
    print("\n  dt   V0=FP-vac      m_eff(V0)      V1=full cumul   m_eff(V1)")
    for dt in range(1, 20):
        print("  %2d  %+11.4e  %7.4f(%.4f)  %+11.4e  %7.4f(%.4f)"
              % (dt, Cc0[dt], em0[dt], ee0[dt], Cc[dt], emm[dt], eme[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    pos = Cc[dts] > 0
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    a1.errorbar(dts[pos], Cc[dts][pos], yerr=err[dts][pos], color="tab:red", marker="o", ms=4, lw=1, capsize=2)
    a1.set_yscale("log")
    a1.set_xlabel(r"$dt$")
    a1.set_ylabel(r"$\langle\sigma\sigma\sigma\sigma\rangle_c(dt)$")
    a1.set_title("connected four-point (log)")
    edt = np.arange(1, twin - 1)
    good = np.isfinite(emm[edt]) & (Cc[edt] > 0) & (Cc[edt + 1] > 0)
    a2.errorbar(edt[good], emm[edt][good], yerr=eme[edt][good], color="tab:red", marker="o", ms=5, lw=1, capsize=2)
    a2.axhline(atm_f2, color="tab:green", ls="-", lw=1.6, alpha=0.85, label=r"0++ $F^2$: $a_t m=%.3f$" % atm_f2)
    a2.axhline(0.67, color="gray", ls="-.", lw=1, alpha=0.5, label=r"$2 a_t m_\sigma\approx0.67$")
    a2.set_ylim(0.0, 1.6)
    a2.set_xlabel(r"$dt$")
    a2.set_ylabel(r"$a_t m_\mathrm{eff}$")
    a2.set_title("connected four-point effmass")
    a2.legend(fontsize=9)
    fig.suptitle("Connected 4-pt cumulant of sigma=PS-1/2  (%s, %d cfg)" % (tag, ncfg))
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_connected_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
