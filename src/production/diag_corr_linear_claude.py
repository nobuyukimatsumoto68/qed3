#!/usr/bin/env python3
# diag_corr_linear_claude.py  [per-diagram CORRELATOR of <PS^2 PS^2>, LINEAR scale, plateau (t-sum) subtracted]
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 \
#         NVDIR=distill_Nv24 CHANNEL=PS CONTACT=0.5 python3 diag_corr_linear_claude.py
#
# The 10 diagrams A..J of <PS^2 PS^2> (diags_pair), translation-averaged, with the large-t plateau
# (t-sum) subtracted per diagram (Cc = C - <C>_{dt>=PLAT_LO}) and config jackknife.  Plots the CORRELATOR
# C_i(dt) (not the effmass) on a LINEAR y-axis.  Same pieces as diag_effmass_claude.py.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

PLAT_LO = int(os.environ.get("PLAT_LO", str(de.PLAT_LO)))
DTMAX = int(os.environ.get("DTMAX", "16"))


def main():
    de.CONTACT = float(os.environ.get("CONTACT", "0.5"))
    tag = dc.ENS.split("nu0")[0]
    channel = os.environ.get("CHANNEL", "PS")
    print("# ENS=%s ncfg=%d  per-diagram CORRELATOR (linear)  CHANNEL=%s CONTACT=%.2f PLAT_LO=%d"
          % (tag, len(dc.KS), channel, de.CONTACT, PLAT_LO))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        leg = tau if channel == "PS" else -taugw
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += de.diags_pair(Phi, leg, s, s + dt)
            D[:, dt] = (acc / ns).real
        allD.append(D)
    allD = np.array(allD)                                  # (ncfg, 10, twin)
    ncfg, ndg, twin = allD.shape

    def jk_corr(C):                                        # C (ncfg, twin) -> plateau-subtracted mean, err
        plat = C[:, PLAT_LO:].mean(1, keepdims=True)
        Cc = C - plat
        n = C.shape[0]
        samp = np.array([np.delete(Cc, i, 0).mean(0) for i in range(n)])
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    corr = {}
    for i in range(10):
        corr[i] = jk_corr(allD[:, i, :])
    Ssum = 2.0 * np.tensordot(dc.W10, allD, axes=(0, 1))
    cT, eT = jk_corr(Ssum)

    print("\n#  dt |  " + "  ".join("%-9s" % de.LABELS[i].split("(")[0] for i in range(10)) + "  TOTAL")
    for dt in range(1, min(DTMAX, twin)):
        row = "  ".join("%9.2e" % corr[i][0][dt] for i in range(10))
        print("#  %2d | %s  %9.2e" % (dt, row, cT[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    panels = list(range(10)) + ["T"]                        # 10 diagrams + TOTAL
    fig, axs = plt.subplots(3, 4, figsize=(14, 9))
    axs = axs.ravel()
    for p, key in enumerate(panels):
        ax = axs[p]
        if key == "T":
            cm, ee, lab, col = cT, eT, "TOTAL (PS.PS)", "black"
        else:
            cm, ee = corr[key]
            lab, col = de.LABELS[key], "tab:blue"
        ax.errorbar(dts, cm[dts], yerr=ee[dts], color=col, marker="o", ms=4, lw=1, capsize=2)
        ax.axhline(0.0, color="gray", lw=0.8, alpha=0.6)
        ax.set_title(lab, fontsize=10)
        ax.set_xlabel(r"$dt$", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[-1].axis("off")
    fig.suptitle("Per-diagram correlator (linear, plateau-subtracted)  <PS^2 PS^2>  %s L1 %d cfg contact=%.2f"
                 % (tag, ncfg, de.CONTACT), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_corr_linear_c%.1f_%s_claude.png" % (de.CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
