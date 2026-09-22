#!/usr/bin/env python3
# fs_diag_corr_v2_claude.py  [CORRECTED FS diagram-by-diagram <sigma_FS^2 sigma_FS^2>, per-loop S/Stilde rule]
# Run:  ENS=... NVDIR=distill_Nv24 python3 fs_diag_corr_v2_claude.py
#
# FIX (NM 2026-09-09): the S/Stilde kernels sum PER CLOSED LOOP (sum first, then multiply loops), NOT per
# diagram; AND always use the 1/2-subtracted (improved) overlap propagator tilde_tau = tau - 1/2 I.
#   - By GW the FS backward-furnished leg collapses: -(1-D_ov^dag) D^{-dag} = D^{-1}, so each loop's
#     (S + Stilde) = 2 x (improved loop).  Hence FS diagram_i = 2^{n_loops[i]} x diags_pair(tilde_tau)[i].
#   - The improved propagator makes every SINGLE-propagator (tadpole) loop vanish: tr[Phi tilde_tau]=0.
#     => diagrams carrying a D_S tadpole (C,D,G,H,I,J) COLLAPSE TO ZERO; only A,B,E (and vacuum F) survive.
# This replaces the old fs_diag_corr_linear (which used leg=-taugw with the wrong FS contact -> spurious
# D_S^FS=-0.96 and a large G).  See fs_sigma2_diagram_note_claude.md.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

PLAT_LO = int(os.environ.get("PLAT_LO", str(de.PLAT_LO)))
DTMAX = int(os.environ.get("DTMAX", "16"))
# closed-loop count per diagram A..J (A=S_S 1 loop; C=D_S*V_S 2; E=C_S^2 2; G=D_S^2 C_S 3; J=D_S^4 4)
NLOOPS = np.array([1, 1, 2, 2, 2, 2, 3, 3, 3, 4])
FSFAC = 2.0 ** NLOOPS


def main():
    de.CONTACT = 0.5                                   # ALWAYS the improved (1/2-subtracted) propagator
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s ncfg=%d  CORRECTED FS <sigma_FS^2 sigma_FS^2> (per-loop S/Stilde, improved prop)" % (tag, len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    tad = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Iv = np.eye(tau.shape[-1])
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        tad.append(np.mean([np.trace(Phi[a] @ (tau[a, a] - 0.5 * Iv)).real for a in range(twin)]))
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += FSFAC * de.diags_pair(Phi, tau, s, s + dt)   # improved forward + per-loop FS factor
            D[:, dt] = (acc / ns).real
        allD.append(D)
    allD = np.array(allD)
    ncfg, ndg, twin = allD.shape
    print("# improved tadpole  <D_S> = %.3e +- %.1e  (== 0 => C,D,G,H,I,J vanish)"
          % (np.mean(tad), np.std(tad) / np.sqrt(ncfg)))

    TSUM_LO = int(os.environ.get("TSUM_LO", "1"))       # t-sum (t-mean) over dt in [TSUM_LO, twin), PER CONFIG
    # PER-CONFIG t-sum subtraction (removes each config's DC/zero-mode; changes disconnected diagrams + errors)
    allD_ts = allD - allD[:, :, TSUM_LO:].mean(axis=2, keepdims=True)
    # per-diagram plateau window [lo,hi) for the SECOND subtraction (identified per diagram from the flat tail)
    DEFPLAT = (twin - 8, twin)
    PLATWIN = {i: DEFPLAT for i in range(10)}

    def dsub_jk(Di, plo, phi):
        # Di = per-config t-sum-subtracted correlator (ncfg, twin).  jackknife + per-diagram plateau subtraction.
        n = Di.shape[0]
        samp = np.array([np.delete(Di, k, 0).mean(0) for k in range(n)])
        samp = samp - samp[:, plo:phi].mean(1, keepdims=True)
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    # inspect the (per-config) t-sum-subtracted MEAN tails to choose per-diagram plateau windows
    print("\n# t-sum-subtracted MEAN tail (pick per-diagram plateau window where flat):")
    print("#  dt |  " + "  ".join("%-9s" % de.LABELS[i].split("(")[0] for i in range(10)))
    tsmean = allD_ts.mean(0)
    for dt in list(range(twin - 10, twin)):
        print("#  %2d | %s" % (dt, "  ".join("%9.2e" % tsmean[i, dt] for i in range(10))))

    # per-diagram DOUBLY-subtracted jackknife samples (t-sum per config -> jk-avg -> per-diagram plateau)
    n = allD_ts.shape[0]
    samp_diag = np.zeros((n, 10, twin))
    for k in range(n):
        Djk = np.delete(allD_ts, k, 0).mean(0)                       # (10, twin)
        for i in range(10):
            plo, phi = PLATWIN[i]
            samp_diag[k, i] = Djk[i] - Djk[i, plo:phi].mean()
    corr = {i: (samp_diag[:, i].mean(0),
                np.sqrt((n - 1) * np.mean((samp_diag[:, i] - samp_diag[:, i].mean(0)) ** 2, 0)))
            for i in range(10)}
    # TOTAL: sum the doubly-subtracted diagrams (W10), then a THIRD plateau subtraction on the total, per sample
    samp_tot = np.einsum('i,kit->kt', dc.W10, samp_diag)
    tlo, thi = DEFPLAT
    samp_tot = samp_tot - samp_tot[:, tlo:thi].mean(1, keepdims=True)
    cT = samp_tot.mean(0)
    eT = np.sqrt((n - 1) * np.mean((samp_tot - cT) ** 2, 0))

    print("\n#  dt |  " + "  ".join("%-9s" % de.LABELS[i].split("(")[0] for i in range(10)) + "  TOTAL")
    for dt in range(1, min(DTMAX, twin)):
        row = "  ".join("%9.2e" % corr[i][0][dt] for i in range(10))
        print("#  %2d | %s  %9.2e" % (dt, row, cT[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    panels = list(range(10)) + ["T"]
    fig, axs = plt.subplots(3, 4, figsize=(14, 9))
    axs = axs.ravel()
    for p, key in enumerate(panels):
        ax = axs[p]
        if key == "T":
            cm, ee, lab, col = cT, eT, "TOTAL (FS.FS)", "black"
        else:
            cm, ee = corr[key]
            lab, col = de.LABELS[key], "tab:purple"
        ax.errorbar(dts, cm[dts], yerr=ee[dts], color=col, marker="o", ms=4, lw=1, capsize=2)
        ax.axhline(0.0, color="gray", lw=0.8, alpha=0.6)
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel(r"$dt$", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[-1].axis("off")
    fig.suptitle("CORRECTED FS per-diagram (per-loop S/Stilde, improved prop): C,D,G,H,I,J vanish  %s L1 %d cfg"
                 % (tag, ncfg), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_diag_corr_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)

    # log |C| version
    figL, axsL = plt.subplots(3, 4, figsize=(14, 9))
    axsL = axsL.ravel()
    for p, key in enumerate(panels):
        ax = axsL[p]
        if key == "T":
            cm, ee, lab, col = cT, eT, "TOTAL (FS.FS)", "black"
        else:
            cm, ee = corr[key]
            lab, col = de.LABELS[key], "tab:purple"
        ax.errorbar(dts, np.abs(cm[dts]), yerr=ee[dts], color=col, marker="o", ms=4, lw=1, capsize=2)
        ax.set_yscale("log")
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel(r"$dt$", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3, which="both")
    axsL[-1].axis("off")
    figL.suptitle("CORRECTED FS per-diagram |correlator| (LOG, double-subtracted)  %s L1 %d cfg"
                  % (tag, ncfg), fontsize=12)
    figL.tight_layout(rect=[0, 0, 1, 0.97])
    outL = "figs/fs_diag_corr_v2_log_%s_claude.png" % tag
    figL.savefig(outL, dpi=130)
    plt.close(figL)
    print("# -> %s" % outL)

    # single-panel log|TOTAL| (triply-subtracted FS.FS correlator)
    figT, axT = plt.subplots(figsize=(8, 5.4))
    axT.errorbar(dts, np.abs(cT[dts]), yerr=eT[dts], color="black", marker="o", ms=5, lw=1.2, capsize=2.5)
    axT.set_yscale("log")
    axT.set_xlabel(r"$dt$")
    axT.set_ylabel(r"$|\langle\sigma_{FS}^2\,\sigma_{FS}^2\rangle_c|$")
    axT.set_title("FS.FS TOTAL |correlator| (LOG, triply-subtracted: t-sum + per-diagram plateau + total plateau)  %s L1 %d cfg"
                  % (tag, ncfg), fontsize=10.5)
    axT.grid(alpha=0.3, which="both")
    figT.tight_layout()
    outT = "figs/fs_diag_corr_v2_totlog_%s_claude.png" % tag
    figT.savefig(outT, dpi=130)
    plt.close(figT)
    print("# -> %s" % outT)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
