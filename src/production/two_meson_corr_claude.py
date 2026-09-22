#!/usr/bin/env python3
# two_meson_corr_claude.py
#
# The VANILLA two-meson correlator <sigma_a sigma_a(t) sigma_a sigma_a(0)>_conn for a in {PS,FS}, from the
# EXACT distillation four-point (10 diagrams A-J), done PROPERLY:
#   - TRANSLATION-AVERAGED over the window source timeslices s (source s, sink s+dt), using the all-to-all
#     perambulator within the window.
#   - vacuum-subtracted by the large-dt plateau (= <sigma sigma>^2, self-consistent in the 4pt normalization).
#   - effective mass with config jackknife.
# NO GEVP, NO basis proliferation -- just the raw two-meson correlator and its effmass.
#
# PS.PS = 2 G10[tau] ; FS.FS = G10[tau] + G10[-tau_gw]  (S_4 + S~_4).  Reuses distill_contract_claude.py.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

PLAT_LO = 24            # plateau region [PLAT_LO, twin) for the vacuum subtraction (connected part dead)


def S4_pair(Phi, leg, si, ti):
    # weighted 10-diagram sum for source si, sink ti (single leg object = one S_4/S~_4 term)
    Ps = Phi[si]
    Pt = Phi[ti]
    tss = leg[si, si]
    ttt = leg[ti, ti]
    tst = leg[si, ti]
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
    diags = np.array([-SS_st, -TS, DSt * VS_st, DSs * VS_ts, CS ** 2,
                      DpSs * DpSt, -DSs * DSt * CS, -DSs ** 2 * DpSt, -DSt ** 2 * DpSs,
                      DSs ** 2 * DSt ** 2])
    return (dc.W10 * diags).sum()


def fourpoint_ta(Phi, leg, twin):
    # translation-averaged G10[leg](dt) over window sources
    G = np.zeros(twin, complex)
    for dt in range(twin):
        ns = twin - dt
        acc = 0.0
        for s in range(ns):
            acc += S4_pair(Phi, leg, s, s + dt)
        G[dt] = acc / ns
    return G


def main():
    print("# ENS=%s  ncfg=%d  (vanilla two-meson correlator)" % (dc.ENS.split("nu0")[0], len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    PS = []
    FS = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        Gt = fourpoint_ta(Phi, tau, twin).real
        # FS leg collapses to tau by GW (S~ D_ov^{-dag} = tau); -taugw was the tau_gw artifact -> FS == PS.
        # See fs_furnishing_derivation_claude.md.  Original (A/B): Gg = fourpoint_ta(Phi, -taugw, twin).real
        Gg = fourpoint_ta(Phi, tau, twin).real
        PS.append(2.0 * Gt)               # PS.PS = 2 G10[tau]
        FS.append(Gt + Gg)                # FS.FS = G10[tau] + G10[-tau']
    PS = np.array(PS)                     # (ncfg, twin)
    FS = np.array(FS)
    ncfg, twin = PS.shape

    def connected(C):
        # C (ncfg,twin) raw 4pt; return per-config connected = C - plateau (plateau per config)
        plat = C[:, PLAT_LO:].mean(1, keepdims=True)
        return C - plat

    def jk_effmass(Craw):
        Cc = connected(Craw)
        n = Cc.shape[0]
        samp = []
        for i in range(n):
            mask = np.ones(n, bool)
            mask[i] = False
            samp.append(Cc[mask].mean(0))
        samp = np.array(samp)                     # (n, twin)
        with np.errstate(all="ignore"):
            em = np.log(samp[:, :-1] / samp[:, 1:])
        emean = em.mean(0)
        eerr = np.sqrt((n - 1) * np.mean((em - emean) ** 2, 0))
        cmean = samp.mean(0)
        return cmean, emean, eerr

    res = {}
    for name, C in (("PS.PS", PS), ("FS.FS", FS)):
        cmean, emean, eerr = jk_effmass(C)
        res[name] = (emean, eerr)
        print("\n=== %s vanilla two-meson correlator (connected, translation-avg) ===" % name)
        print("  dt    C_conn        m_eff(err)")
        for dt in range(1, 20):
            print("  %2d   %11.4e   %7.4f(%.4f)" % (dt, cmean[dt], emean[dt], eerr[dt]))

    # ---- effective-mass plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tag = dc.ENS.split("nu0")[0]
    dts = np.arange(1, 18)
    fig, ax = plt.subplots(figsize=(7, 5))
    # two-meson effmasses: distinct color AND marker (color-blind safe)
    emp, eep = res["PS.PS"]
    emf, eef = res["FS.FS"]
    ax.errorbar(dts, emp[dts], yerr=eep[dts], color="tab:red", marker="o", ms=5, lw=1, capsize=2,
                label=r"PS$\cdot$PS  two-meson")
    ax.errorbar(dts + 0.1, emf[dts], yerr=eef[dts], color="tab:blue", marker="s", ms=5, lw=1, capsize=2,
                label=r"FS$\cdot$FS  two-meson")
    # single-meson C_S reference masses (from the C_S plateaus) and the 2 m_sigma threshold
    ax.axhline(0.35, color="tab:red", ls="--", lw=1, alpha=0.7, label=r"single $\sigma_{PS}$ ($m\approx0.35$)")
    ax.axhline(0.32, color="tab:blue", ls=":", lw=1, alpha=0.7, label=r"single $\sigma_{FS}$ ($m\approx0.32$)")
    ax.axhline(0.67, color="k", ls="-.", lw=1, alpha=0.6, label=r"$2 m_\sigma\approx0.67$ (two-particle)")
    ax.set_xlabel(r"$dt$  (a_t units)")
    ax.set_ylabel(r"$a_t\, m_\mathrm{eff}$")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Vanilla two-meson effective mass (%s, L1, 400 cfg)" % tag)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_effmass_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
