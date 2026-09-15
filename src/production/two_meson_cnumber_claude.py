#!/usr/bin/env python3
# two_meson_cnumber_claude.py
# CORRECT sigma = PS - 1/2 four-point: the -1/2 is a C-NUMBER (contact VEV) subtraction, NOT a kernel
# modification.  So it removes ONLY the standalone tadpole D_S (self-contraction):  D_S -> D_S - c0,
# c0 = (1/2) Tr[Phi].  The equal-time tau(s,s) INSIDE the connected loops (D'_S, V_S, S_S) stays RAW
# (it links two different sigma's, not a self-contraction).  Contrast with the old (wrong) tau(s,s)-1/2 I
# in every leg (`diag_effmass_claude.py`).
#
# Key consequence: the F diagram D'_S(s) D'_S(t) is now RAW and NONZERO (the old code killed it to noise).
# Its CONNECTED part <dD'_S(s) dD'_S(t)>_c is the disconnected-same-timeslice local sigma^2 operator =
# the fermionic 0++ candidate (compare a_t m_F2 = 0.616).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de       # OLD (wrong) diags_pair for comparison

LABELS = de.LABELS
at = 0.2
atm_f2 = at * 3.08


def diags_pair_cnumber(Phi, leg, si, ti):
    Ps = Phi[si]
    Pt = Phi[ti]
    tss = leg[si, si]                          # RAW equal-time
    ttt = leg[ti, ti]                          # RAW
    tst = leg[si, ti]
    tts = leg[ti, si]
    c0s = 0.5 * np.trace(Ps)                    # tadpole VEV at s  = (1/2) Tr[Phi]
    c0t = 0.5 * np.trace(Pt)
    DSs = np.trace(Ps @ tss) - c0s             # normal-ordered tadpole (self-contraction only)
    DSt = np.trace(Pt @ ttt) - c0t
    DpSs = np.trace(Ps @ tss @ Ps @ tss)       # RAW connected equal-time loop
    DpSt = np.trace(Pt @ ttt @ Pt @ ttt)
    M = Ps @ tst @ Pt @ tts
    CS = np.trace(M)
    TS = np.trace(M @ M)
    VS_st = np.trace(Ps @ tss @ M)             # RAW tss inside V_S
    VS_ts = np.trace(Pt @ ttt @ Pt @ tts @ Ps @ tst)
    SS_st = np.trace(Ps @ tss @ Ps @ tst @ Pt @ ttt @ Pt @ tts)   # RAW inside S_S
    return np.array([-SS_st, -TS, DSt * VS_st, DSs * VS_ts, CS ** 2,
                     DpSs * DpSt, -DSs * DSt * CS, -DSs ** 2 * DpSt, -DSt ** 2 * DpSs,
                     DSs ** 2 * DSt ** 2])


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  sigma=PS-1/2 (c-number tadpole-only) four-point" % (tag, len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []            # corrected per-diagram (ncfg,10,twin)
    DpS = []             # RAW D'_S(s) per config,timeslice for the F-connected 0++ test
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += diags_pair_cnumber(Phi, tau, s, s + dt)
            D[:, dt] = (acc / ns).real
        allD.append(D)
        dps = np.array([np.trace(Phi[s] @ tau[s, s] @ Phi[s] @ tau[s, s]).real for s in range(twin)])
        DpS.append(dps)
    allD = np.array(allD)
    DpS = np.array(DpS)                          # (ncfg, twin)
    ncfg, _, twin = allD.shape

    # ---- per-diagram effmass (plateau-subtracted, config jk) ----
    def jk_effmass(C):
        Cc = C - C[:, 24:].mean(1, keepdims=True)
        s = np.array([np.delete(Cc, i, 0).mean(0) for i in range(ncfg)])
        with np.errstate(all="ignore"):
            em = np.log(s[:, :-1] / s[:, 1:])
        return em.mean(0), np.sqrt((ncfg - 1) * np.mean((em - em.mean(0)) ** 2, 0)), s.mean(0)

    print("\n  diagram        m(dt4)   m(dt8)   m(dt12)   |C(dt4)|   (c-number, RAW connected)")
    for i in range(10):
        em, ee, cm = jk_effmass(allD[:, i, :])
        print("  %-14s  %6.3f   %6.3f   %6.3f   %.3e" % (LABELS[i], em[4], em[8], em[12], abs(cm[4])))

    # ---- F-connected 0++ test:  <dD'_S(s) dD'_S(s+dt)>_c ,  dD'_S = D'_S - <D'_S> ----
    print("\n  <D'_S>_raw = %.5f   (was -0.00985 with the WRONG tau-1/2 I)" % DpS.mean())
    GF = np.zeros((ncfg, twin))
    for dt in range(twin):
        ns = twin - dt
        acc = np.zeros(ncfg)
        for s in range(ns):
            acc += DpS[:, s] * DpS[:, s + dt]
        GF[:, dt] = acc / ns
    # connected + jackknife effmass
    samp = np.zeros((ncfg, twin))
    for i in range(ncfg):
        m = np.ones(ncfg, bool)
        m[i] = False
        samp[i] = GF[m].mean(0) - DpS[m].mean() ** 2
    Cc = samp.mean(0)
    err = np.sqrt((ncfg - 1) * np.mean((samp - Cc) ** 2, 0))
    with np.errstate(all="ignore"):
        em = np.log(samp[:, :-1] / samp[:, 1:])
    emm = em.mean(0)
    eme = np.sqrt((ncfg - 1) * np.mean((em - emm) ** 2, 0))
    print("\n  F-connected  <dD'_S dD'_S>_c   (fermionic 0++ candidate; 0++ a_t m=0.616)")
    print("  dt   C_conn        |C/err|   m_eff(err)")
    for dt in range(1, 16):
        print("  %2d  %+11.4e  %7.2f   %7.4f(%.4f)" % (dt, Cc[dt], abs(Cc[dt]) / err[dt], emm[dt], eme[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    edt = np.arange(1, twin - 1)
    good = np.isfinite(emm[edt]) & (Cc[edt] > 0) & (Cc[edt + 1] > 0)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.errorbar(edt[good], emm[edt][good], yerr=eme[edt][good], color="tab:red", marker="o", ms=5, lw=1,
                capsize=2, label=r"F-connected $\langle\delta D'_S\,\delta D'_S\rangle_c$ (local $\sigma^2$)")
    ax.axhline(atm_f2, color="tab:green", ls="-", lw=1.6, alpha=0.85, label=r"0++ $F^2$: $a_t m=%.3f$" % atm_f2)
    ax.axhline(0.67, color="gray", ls="-.", lw=1, alpha=0.5, label=r"$2 a_t m_\sigma\approx0.67$")
    ax.set_ylim(0.0, 1.4)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("Disconnected same-timeslice (F, RAW) connected corr -- 0++ test (%s, %d cfg)" % (tag, ncfg))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_cnumber_Fconn_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
