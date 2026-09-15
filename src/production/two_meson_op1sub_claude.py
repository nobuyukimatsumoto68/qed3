#!/usr/bin/env python3
# two_meson_op1sub_claude.py
# Two-meson PS.PS summed effmass with the FULL data-measured equal-time one-point subtracted, instead of
# just the analytic contact 1/2.  Two-pass:
#   pass 1: M(a) = <tau(a,a)>_cfg  (Nv x Nv matrix, per window-time a) -- the measured one-point.
#   pass 2: rebuild the 10 diagrams with tau(s,s) -> tau(s,s) - M(s)  (and same at t), sum the kept
#           diagrams (drop F,H,I noise), per-jackknife plateau subtraction, effmass + 0++ line.
# Compared side by side against the scalar-contact (1/2 I) subtraction.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

DROP = [5, 7, 8]                                       # F, H, I noise
KEEP = [i for i in range(10) if i not in DROP]
PLAT_LO = 24
at = 0.2
atm_f2 = at * 3.08                                     # 0++ line in a_t m units


def diags10(Phi_s, Phi_t, tss, ttt, tst, tts):
    DSs = np.trace(Phi_s @ tss)
    DSt = np.trace(Phi_t @ ttt)
    DpSs = np.trace(Phi_s @ tss @ Phi_s @ tss)
    DpSt = np.trace(Phi_t @ ttt @ Phi_t @ ttt)
    M = Phi_s @ tst @ Phi_t @ tts
    CS = np.trace(M)
    TS = np.trace(M @ M)
    VS_st = np.trace(Phi_s @ tss @ M)
    VS_ts = np.trace(Phi_t @ ttt @ Phi_t @ tts @ Phi_s @ tst)
    SS_st = np.trace(Phi_s @ tss @ Phi_s @ tst @ Phi_t @ ttt @ Phi_t @ tts)
    return np.array([-SS_st, -TS, DSt * VS_st, DSs * VS_ts, CS ** 2,
                     DpSs * DpSt, -DSs * DSt * CS, -DSs ** 2 * DpSt, -DSt ** 2 * DpSs,
                     DSs ** 2 * DSt ** 2])


def summed_corr(sub_mode, M=None):
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    wkeep = 2.0 * dc.W10[KEEP]
    C = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Iv = np.eye(tau.shape[-1])
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        c = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                t = s + dt
                tss = tau[s, s] - (0.5 * Iv if sub_mode == "contact" else M[s])
                ttt = tau[t, t] - (0.5 * Iv if sub_mode == "contact" else M[t])
                acc += diags10(Phi[s], Phi[t], tss, ttt, tau[s, t], tau[t, s])
            acc = (acc / ns).real
            c[dt] = float(wkeep @ acc[KEEP])
        C.append(c)
    return np.array(C)


def measure_M():
    Msum = None
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        d = np.array([tau[a, a] for a in range(twin)])       # (twin, Nv, Nv)
        Msum = d if Msum is None else Msum + d
    return Msum / len(dc.KS)


def effmass(C):
    ncfg, twin = C.shape
    samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])
    samp = samp - samp[:, PLAT_LO:].mean(1, keepdims=True)
    with np.errstate(all="ignore"):
        em = np.log(samp[:, :-1] / samp[:, 1:])
    return samp.mean(0), em.mean(0), np.sqrt((ncfg - 1) * np.mean((em - em.mean(0)) ** 2, 0))


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  two-meson: contact vs measured one-point subtraction" % (tag, len(dc.KS)))
    M = measure_M()
    Cc = summed_corr("contact")
    Cm = summed_corr("matrix", M)
    cm_c, em_c, ee_c = effmass(Cc)
    cm_m, em_m, ee_m = effmass(Cm)
    twin = Cc.shape[1]
    print("\n  dt   m_eff(contact)     m_eff(<tau> one-point)")
    for dt in range(1, 20):
        print("  %2d   %7.4f(%.4f)     %7.4f(%.4f)" % (dt, em_c[dt], ee_c[dt], em_m[dt], ee_m[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    edt = np.arange(1, twin - 1)
    gc = np.isfinite(em_c[edt]) & (cm_c[edt] > 0) & (cm_c[edt + 1] > 0)
    gm = np.isfinite(em_m[edt]) & (cm_m[edt] > 0) & (cm_m[edt + 1] > 0)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.errorbar(edt[gc], em_c[edt][gc], yerr=ee_c[edt][gc], color="tab:blue", marker="s", ms=4, lw=1,
                capsize=2, alpha=0.7, label="contact (1/2 I) subtraction")
    ax.errorbar(edt[gm] + 0.12, em_m[edt][gm], yerr=ee_m[edt][gm], color="tab:red", marker="o", ms=5, lw=1,
                capsize=2, label=r"measured one-point $\langle\tau(s,s)\rangle$ subtraction")
    ax.axhline(atm_f2, color="tab:green", ls="-", lw=1.6, alpha=0.85, label=r"0++ $F^2$: $a_t m=%.3f$" % atm_f2)
    ax.axhline(0.67, color="gray", ls="-.", lw=1, alpha=0.5, label=r"$2 a_t m_\sigma\approx0.67$")
    ax.set_ylim(0.0, 1.4)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    ax.set_title("Two-meson effmass: contact vs measured one-point subtraction (%s, %d cfg)" % (tag, len(dc.KS)))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_op1sub_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
