#!/usr/bin/env python3
# fs_diag_corr_linear_claude.py  [per-diagram CORRELATOR of <FS^2 FS^2>, LINEAR, plateau (t-sum) subtracted]
# Run:  ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 \
#         NVDIR=distill_Nv24 python3 fs_diag_corr_linear_claude.py
#
# FS channel (see fs_sigma2_diagram_note_claude.md):  FS.FS = sum_i W10[i] ( G10[tau]_i + G10[-tau']_i ),
# NOT 2 G10[-tau'].  Two parts per diagram:
#   S-part   : leg = tau      , equal-time contact tt(s,s)=tau(s,s) - 1/2 I         (L1: O=V^d V=I)
#   Stilde   : leg = -tau'    , equal-time contact tt(s,s)=-1/2 ( tau'(s,s)+tau(s,s) )
#              (furnished factor (1-D_ov^d) multiplies the 1/2 overlap contact: contact(tau')=1/2(tau'-tau))
#   off-diagonal legs carry NO contact: tau(s,t) for S ; -tau'(s,t) for Stilde.
# Per-diagram FS_i = S_i + Stilde_i ; TOTAL = sum_i W10[i] (S_i + Stilde_i).  Plateau (t-sum) subtracted,
# config jackknife, LINEAR per-diagram panels.  tau'=tau_gw (furnished perambulator).

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


def diags_from_legs(Ps, Pt, tss, ttt, tst, tts):
    # the 10 diagrams from explicit equal-time (tss,ttt) and off-diagonal (tst,tts) legs
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
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s ncfg=%d  per-diagram CORRELATOR (linear)  CHANNEL=FS  PLAT_LO=%d" % (tag, len(dc.KS), PLAT_LO))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    tad_S = []
    tad_St = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Iv = np.eye(tau.shape[-1])
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        # equal-time contact-subtracted diagonals
        ttS = [tau[a, a] - 0.5 * Iv for a in range(twin)]                     # S-part  (leg tau)
        ttSt = [-0.5 * (taugw[a, a] + tau[a, a]) for a in range(twin)]        # Stilde  (leg -tau')
        # tadpole diagnostic <Tr[Phi tt]>
        tad_S.append(np.mean([np.trace(Phi[a] @ ttS[a]).real for a in range(twin)]))
        tad_St.append(np.mean([np.trace(Phi[a] @ ttSt[a]).real for a in range(twin)]))
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                t = s + dt
                # S-part: off-diagonals tau(s,t), tau(t,s)
                sS = diags_from_legs(Phi[s], Phi[t], ttS[s], ttS[t], tau[s, t], tau[t, s])
                # Stilde-part COLLAPSES to the S-part by GW (S~ D_ov^{-dag} = tau): use ttS/tau legs, NOT
                # ttSt/-taugw (the tau_gw artifact).  See fs_furnishing_derivation_claude.md -> FS == PS.
                # Original (A/B): sSt = diags_from_legs(Phi[s], Phi[t], ttSt[s], ttSt[t], -taugw[s, t], -taugw[t, s])
                sSt = diags_from_legs(Phi[s], Phi[t], ttS[s], ttS[t], tau[s, t], tau[t, s])
                acc += sS + sSt
            D[:, dt] = (acc / ns).real
        allD.append(D)
    allD = np.array(allD)                                  # (ncfg, 10, twin)  = G10[tau]+G10[-tau'] per diagram
    ncfg, ndg, twin = allD.shape
    print("# FS tadpole <Tr[Phi tt]>  S-part(tau) = %.3e +- %.1e ;  Stilde(-tau') = %.3e +- %.1e"
          % (np.mean(tad_S), np.std(tad_S) / np.sqrt(ncfg),
             np.mean(tad_St), np.std(tad_St) / np.sqrt(ncfg)))

    def jk_corr(C):
        plat = C[:, PLAT_LO:].mean(1, keepdims=True)
        Cc = C - plat
        n = C.shape[0]
        samp = np.array([np.delete(Cc, i, 0).mean(0) for i in range(n)])
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    corr = {i: jk_corr(allD[:, i, :]) for i in range(10)}
    Ssum = np.tensordot(dc.W10, allD, axes=(0, 1))         # sum_i W10 (S_i + Stilde_i)  (NO factor 2 for FS)
    cT, eT = jk_corr(Ssum)

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
    fig.suptitle("Per-diagram correlator (linear, plateau-sub)  <FS^2 FS^2>=G10[tau]+G10[-tau']  %s L1 %d cfg"
                 % (tag, ncfg), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_diag_corr_linear_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
