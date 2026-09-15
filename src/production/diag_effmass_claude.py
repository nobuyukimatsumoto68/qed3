#!/usr/bin/env python3
# diag_effmass_claude.py
# Per-diagram effective mass of the 10 two-meson diagrams (PS channel, legs=tau), translation-averaged
# and plateau-subtracted, with config jackknife.  Diagnostic: which diagram carries which mass?

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

LABELS = ["A(-S_S)", "B(-T_S)", "C(D_S V_S)", "D(D_S V_S)", "E(C_S^2)",
          "F(D'D')", "G(-D_S^2 C_S)", "H(-D_S^2 D')", "I(-D_S^2 D')", "J(D_S^4)"]
PLAT_LO = 24
CONTACT = float(os.environ.get("CONTACT", "0.5"))   # GW <psibar psi> contact = 1/2; env CONTACT=0 to disable


def diags_pair(Phi, leg, si, ti):
    Iv = np.eye(leg.shape[-1])
    Ps = Phi[si]
    Pt = Phi[ti]
    tss = leg[si, si] - CONTACT * Iv                 # normal-order: subtract the contact from equal-time loops
    ttt = leg[ti, ti] - CONTACT * Iv
    tst = leg[si, ti]                                # off-diagonal legs: no contact
    tts = leg[ti, si]
    DSs = np.trace(Ps @ tss)
    DSt = np.trace(Pt @ ttt)
    DpSs = np.trace(Ps @ tss @ Ps @ tss)
    DpSt = np.trace(Pt @ ttt @ Pt @ ttt)
    M = Ps @ tst @ Pt @ tts
    CS = np.trace(M)
    TS = np.trace(M @ M)
    VS_st = np.trace(Ps @ tss @ M)
    VS_ts = np.trace(Pt @ ttt @ Pt @ tts @ Ps @ tst)
    SS_st = np.trace(Ps @ tss @ Ps @ tst @ Pt @ ttt @ Pt @ tts)
    return np.array([-SS_st, -TS, DSt * VS_st, DSs * VS_ts, CS ** 2,
                     DpSs * DpSt, -DSs * DSt * CS, -DSs ** 2 * DpSt, -DSt ** 2 * DpSs,
                     DSs ** 2 * DSt ** 2])


def main():
    print("# ENS=%s  ncfg=%d  per-diagram effmass" % (dc.ENS.split("nu0")[0], len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    channel = os.environ.get("CHANNEL", "PS")     # PS -> legs=tau ; FS -> legs=-taugw
    print("# CHANNEL=%s" % channel)
    allD = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        leg = tau if channel == "PS" else -taugw
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += diags_pair(Phi, leg, s, s + dt)
            D[:, dt] = (acc / ns).real
        allD.append(D)
    allD = np.array(allD)                          # (ncfg, 10, twin)
    ncfg, ndg, twin = allD.shape

    def jk_effmass(C):                             # C (ncfg, twin)
        plat = C[:, PLAT_LO:].mean(1, keepdims=True)
        Cc = C - plat
        n = C.shape[0]
        samp = []
        for i in range(n):
            m = np.ones(n, bool)
            m[i] = False
            samp.append(Cc[m].mean(0))
        samp = np.array(samp)
        with np.errstate(all="ignore"):
            em = np.log(samp[:, :-1] / samp[:, 1:])
        return em.mean(0), np.sqrt((n - 1) * np.mean((em - em.mean(0)) ** 2, 0)), samp.mean(0)

    print("# CONTACT (normal-order) subtraction = %.3f  (0 = raw, 0.5 = GW contact)" % CONTACT)
    ems = {}
    print("\n  diagram        m_eff(dt=4)   m_eff(dt=8)   m_eff(dt=12)   |C_conn(dt=4)|")
    for i in range(10):
        em, ee, cm = jk_effmass(allD[:, i, :])
        ems[i] = (em, ee)
        print("  %-14s  %6.3f(%.3f)  %6.3f(%.3f)  %6.3f(%.3f)   %.3e"
              % (LABELS[i], em[4], ee[4], em[8], ee[8], em[12], ee[12], abs(cm[4])))
    # TOTAL two-meson correlator = weighted sum (PS.PS = 2 * sum)
    Ssum = 2.0 * np.tensordot(dc.W10, allD, axes=(0, 1))     # (ncfg, twin)
    em, ee, cm = jk_effmass(Ssum)
    print("\n  [TOTAL two-meson]  dt   C_conn        m_eff(err)")
    for dt in range(1, 16):
        print("                    %2d   %11.4e   %7.4f(%.4f)" % (dt, cm[dt], em[dt], ee[dt]))

    # plot the diagrams with signal (A,B,C,D,E,G); F,H,I,J are flat/vacuum
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tag = dc.ENS.split("nu0")[0]
    sig = [0, 1, 2, 4, 6]      # A, B, C, E, G  (D==C; H,I,J flat)
    styles = {0: ("tab:red", "o"), 1: ("tab:blue", "s"), 2: ("tab:green", "^"),
              4: ("tab:purple", "D"), 6: ("tab:orange", "v")}
    dts = np.arange(1, 16)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for i in sig:
        em, ee = ems[i]
        col, mk = styles[i]
        ax.errorbar(dts + 0.03 * i, em[dts], yerr=ee[dts], color=col, marker=mk, ms=4, lw=1, capsize=2,
                    label=LABELS[i])
    # TOTAL two-meson effmass (thick black)
    emT, eeT, _ = jk_effmass(Ssum)
    ax.errorbar(dts, emT[dts], yerr=eeT[dts], color="black", marker="*", ms=8, lw=1.8, capsize=2,
                label="TOTAL (PS.PS)")
    ax.axhline(0.35, color="gray", ls="--", lw=1, alpha=0.6, label=r"$a_t m_\sigma\approx0.35$")
    ax.axhline(0.67, color="gray", ls="-.", lw=1, alpha=0.6, label=r"$2 a_t m_\sigma\approx0.67$")
    ax.axhline(0.2 * 3.08, color="tab:green", ls="-", lw=1.6, alpha=0.85, label=r"0++ $F^2$: $a_t m=0.616$")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    ax.set_ylim(0.0, 1.2)
    ax.set_title("Per-diagram effmass (PS legs, contact=%.2f, %s, L1, %d cfg)" % (CONTACT, tag, ncfg))
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_effmass_c%.1f_%s_claude.png" % (CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
