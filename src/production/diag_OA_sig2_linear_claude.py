#!/usr/bin/env python3
# diag_OA_sig2_linear_claude.py  [per-diagram CORRELATOR of <O_A PS^2>, LINEAR, plateau (t-sum) subtracted]
# Run:  ENS=... NVDIR=distill_Nv24 CHANNEL=PS CONTACT=0.5 python3 diag_OA_sig2_linear_claude.py
#
# <O_A(t) sigma^2(s)>, O_A = psibar Phi tilde_tau psi (vertex PA=Phi tt), sigma = psibar Phi psi (vertex Phi).
# 6 fermion fields -> 3 propagators; the 6 Wick pairings group into 4 diagram types (multiplicities):
#   Tri   (x2): triangle, connected coupling   = -Tr[PA(t) tau(t,s) Phi(s) tt(s,s) Phi(s) tau(s,t)]
#   Semi  (x2): O_A-sigma 2-loop * sigma tadpole = +Tr[PA(t) tau(t,s) Phi(s) tau(s,t)] * Tr[Phi(s) tt(s,s)]
#   OA_Dp (x1): O_A tadpole * sigma^2 self-loop  = +Tr[PA(t) tt(t,t)] * Tr[Phi(s) tt(s,s) Phi(s) tt(s,s)]
#   Disc  (x1): O_A tadpole * sigma tadpole^2    = -Tr[PA(t) tt(t,t)] * Tr[Phi(s) tt(s,s)]^2
# (Tri validated: 2*Tri == connected C2A of gevp_OA_sig22_free.)  Sign = (-1)^(#loops).
# Plateau (t-sum) subtracted per diagram, config jackknife.  Linear per-diagram panels.

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
LAB = ["Tri (x2)", "Semi (x2)", "OA_Dp (x1)", "Disc (x1)"]
MULT = np.array([2.0, 2.0, 1.0, 1.0])


def main():
    de.CONTACT = float(os.environ.get("CONTACT", "0.5"))
    tag = dc.ENS.split("nu0")[0]
    channel = os.environ.get("CHANNEL", "PS")
    print("# ENS=%s ncfg=%d  per-diagram <O_A PS^2> (linear, plateau-sub)  CONTACT=%.2f" % (tag, len(dc.KS), de.CONTACT))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        leg = tau if channel == "PS" else -taugw
        tt = [leg[a, a] - de.CONTACT * np.eye(leg.shape[-1]) for a in range(twin)]
        PA = [Phi[a] @ tt[a] for a in range(twin)]
        D = np.zeros((4, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(4)
            for s in range(ns):
                t = s + dt
                Ph = Phi[s]
                tss = tt[s]
                oaloop = np.trace(PA[t] @ leg[t, s] @ Ph @ leg[s, t]).real
                dS = np.trace(Ph @ tss).real
                dpS = np.trace(Ph @ tss @ Ph @ tss).real
                dA = np.trace(PA[t] @ tt[t]).real
                tri = -np.trace(PA[t] @ leg[t, s] @ Ph @ tss @ Ph @ leg[s, t]).real
                semi = oaloop * dS
                oa_dp = dA * dpS
                disc = -dA * dS * dS
                acc += np.array([tri, semi, oa_dp, disc])
            D[:, dt] = acc / ns
        allD.append(D)
    allD = np.array(allD)                                  # (ncfg,4,twin)
    ncfg, _, twin = allD.shape

    def jk_corr(C):
        plat = C[:, PLAT_LO:].mean(1, keepdims=True)
        Cc = C - plat
        n = C.shape[0]
        samp = np.array([np.delete(Cc, i, 0).mean(0) for i in range(n)])
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    corr = [jk_corr(allD[:, i, :]) for i in range(4)]
    tot = (MULT[None, :, None] * allD).sum(1)              # weighted total
    cT, eT = jk_corr(tot)

    print("\n#  dt |  " + "  ".join("%-11s" % LAB[i].split()[0] for i in range(4)) + "  TOTAL")
    for dt in range(1, DTMAX):
        row = "  ".join("%11.3e" % corr[i][0][dt] for i in range(4))
        print("#  %2d | %s  %11.3e" % (dt, row, cT[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    fig, axs = plt.subplots(2, 3, figsize=(13, 7))
    axs = axs.ravel()
    for i in range(4):
        cm, ee = corr[i]
        axs[i].errorbar(dts, cm[dts], yerr=ee[dts], color="tab:blue", marker="o", ms=4, lw=1, capsize=2)
        axs[i].axhline(0.0, color="gray", lw=0.8, alpha=0.6)
        axs[i].set_title(LAB[i], fontsize=10)
        axs[i].set_xlabel(r"$dt$", fontsize=8)
        axs[i].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[4].errorbar(dts, cT[dts], yerr=eT[dts], color="black", marker="*", ms=7, lw=1.4, capsize=2)
    axs[4].axhline(0.0, color="gray", lw=0.8, alpha=0.6)
    axs[4].set_title("TOTAL (2Tri+2Semi+OA_Dp+Disc)", fontsize=10)
    axs[4].set_xlabel(r"$dt$", fontsize=8)
    axs[4].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[5].axis("off")
    fig.suptitle("Per-diagram <O_A PS^2> (linear, plateau-sub)  %s L1 %d cfg contact=%.2f"
                 % (tag, ncfg, de.CONTACT), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_OA_sig2_linear_c%.1f_%s_claude.png" % (de.CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
