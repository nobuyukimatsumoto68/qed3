#!/usr/bin/env python3
# f2_sigma2_cross_v2_claude.py
#   CORRECTED F^2 -- sigma^2 CROSS  <F^2(s+dt) sigma^2_00(s)>_c , distillation, Nf2 gsq1.0 L1.
#   Redo of f2_sigma2_cross_diag_claude.py with the TRIPLE SUBTRACTION (correct errors), per
#   fs_diag_corr_v2_claude.py / the corrected sigma machinery.  PS == FS here EXACTLY: F^2 is
#   gluonic so a single sigma^2(s) contracts only into
#       disc  =  D_S(s)^2        (two length-1 tadpole loops -> vanish under improvement -> DEAD)
#       ext   =  D'_S(s)         (one length-2 loop -> EVEN -> furnishing = +1 -> no FS/PS split)
#   with  D_S = Tr[Phi tilde_tau] ,  D'_S = Tr[Phi tilde_tau Phi tilde_tau] ,
#         tilde_tau = tau(s,s) - CONTACT*I  (GW contact 1/2 ; CONTACT=0.5).  sigma^2 total = 2(D_S^2 + D'_S).
#
#   TRIPLE SUBTRACTION (per jackknife, correct errors):
#     (1) per-config t-sum: C_k(dt) - mean_dt C_k(dt)  (removes each config's vacuum <F^2><X> level)
#     (2) per-diagram plateau: subtract the large-dt (vacuum) tail per jk sample per diagram
#     (3) total plateau: sum 2*(disc+ext), reidentify the tail, subtract, jackknife
#
#   Run: ENS=Nf2_gsq1.000000...L1_hb1.000000 NVDIR=distill_Nv24 python3 f2_sigma2_cross_v2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import h5py
import distill_contract_claude as dc

Nt = 128
OP_F = int(os.environ.get("OP_F", "0"))
CONTACT = float(os.environ.get("CONTACT", "0.5"))
DTMAX = int(os.environ.get("DTMAX", "48"))          # plotted range
DTCALC = int(os.environ.get("DTCALC", "64"))        # computed range = Nt/2 (F^2 global -> full reach; vacuum tail)
TSUM_LO = int(os.environ.get("TSUM_LO", "1"))       # per-config t-sum (t-mean) over dt in [TSUM_LO, DTCALC)
PLAT_LO = int(os.environ.get("PLAT_LO", str(DTCALC - 8)))   # plateau (vacuum tail) window start
LABELS = ["disc  $D_S^2$", "ext  $D'_S$", "TOTAL  $2(D_S^2+D'_S)$"]


def gluedir():
    return "data_" + dc.ENS


def load_of(k):
    with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k), "r") as g:
        return np.array(g["O"])[OP_F].astype(float)     # (Nt,)


def per_config(k, w00):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    tt = [tau[a, a] - CONTACT * np.eye(tau.shape[-1]) for a in range(twin)]
    DS = np.array([np.trace(Phi[a] @ tt[a]).real for a in range(twin)])
    DpS = np.array([np.trace(Phi[a] @ tt[a] @ Phi[a] @ tt[a]).real for a in range(twin)])
    of = load_of(k)                                       # (Nt,)
    Cd = np.zeros((2, DTCALC))                             # [disc, ext]
    for dt in range(DTCALC):
        # clean t-fold: C(dt)=C(-dt) for two Hermitian 0++ ops; the backward F^2 at s-dt is FREE (full-Nt).
        idxf = (tsrc0 + np.arange(twin) + dt) % Nt
        idxb = (tsrc0 + np.arange(twin) - dt) % Nt
        ofd = 0.5 * (of[idxf] + of[idxb])                 # F^2 forward+backward around each sigma^2 source
        Cd[0, dt] = np.mean(ofd * (DS ** 2))              # disc  D_S^2
        Cd[1, dt] = np.mean(ofd * DpS)                    # ext   D'_S
    return Cd


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    print("# ENS=%s  matched cfg=%d/%d  OP_F=%d CONTACT=%.2f  (PS==FS: all loops even/tadpole)"
          % (tag, len(ks), len(dc.KS), OP_F, CONTACT))
    allC = np.array([per_config(k, w00) for k in ks])     # (ncfg, 2, DTCALC)
    ncfg = allC.shape[0]

    # (1) per-config t-sum subtraction (removes each config's vacuum <F^2><X> DC level; correct errors)
    allC_ts = allC - allC[:, :, TSUM_LO:].mean(axis=2, keepdims=True)

    DEFPLAT = (PLAT_LO, DTCALC)
    PLATWIN = {0: DEFPLAT, 1: DEFPLAT}

    # inspect the (per-config) t-sum-subtracted MEAN tail to confirm the plateau window is flat
    tsmean = allC_ts.mean(0)
    print("\n# t-sum-subtracted MEAN tail (dt : disc  ext):")
    for dt in range(DTCALC - 10, DTCALC):
        print("#  %2d | %11.3e %11.3e" % (dt, tsmean[0, dt], tsmean[1, dt]))

    # (2) per-diagram doubly-subtracted jackknife (t-sum per config -> jk-avg -> per-diagram plateau)
    n = ncfg
    samp_diag = np.zeros((n, 2, DTCALC))
    for k in range(n):
        Djk = np.delete(allC_ts, k, 0).mean(0)            # (2, DTCALC)
        for i in range(2):
            plo, phi = PLATWIN[i]
            samp_diag[k, i] = Djk[i] - Djk[i, plo:phi].mean()
    corr = {i: (samp_diag[:, i].mean(0),
                np.sqrt((n - 1) * np.mean((samp_diag[:, i] - samp_diag[:, i].mean(0)) ** 2, 0)))
            for i in range(2)}

    # (3) total = 2*(disc+ext), third plateau subtraction per sample
    samp_tot = 2.0 * (samp_diag[:, 0] + samp_diag[:, 1])
    tlo, thi = DEFPLAT
    samp_tot = samp_tot - samp_tot[:, tlo:thi].mean(1, keepdims=True)
    cT = samp_tot.mean(0)
    eT = np.sqrt((n - 1) * np.mean((samp_tot - cT) ** 2, 0))

    print("\n#  dt |   disc(D_S^2)        ext(D'_S)         TOTAL 2(disc+ext)     ext S/N")
    for dt in range(0, DTMAX):
        sn = corr[1][0][dt] / corr[1][1][dt] if corr[1][1][dt] > 0 else 0.0
        print("#  %2d | %11.3e(%.1e) %11.3e(%.1e) %11.3e(%.1e)  %6.2f"
              % (dt, corr[0][0][dt], corr[0][1][dt], corr[1][0][dt], corr[1][1][dt],
                 cT[dt], eT[dt], sn))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(0, DTMAX)
    panels = [corr[0], corr[1], (cT, eT)]

    # linear
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
    for ax, cc, lab, col, mk in zip(axs, panels, LABELS,
                                    ["tab:red", "tab:green", "black"], ["o", "^", "*"]):
        ax.errorbar(dts, cc[0][dts], yerr=cc[1][dts], color=col, marker=mk, ms=5, lw=1, capsize=2)
        ax.axhline(0.0, color="gray", lw=0.8, alpha=0.6)
        ax.set_title(lab, fontsize=11)
        ax.set_xlabel(r"$dt$", fontsize=9)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[0].set_ylabel(r"$\langle F^2(s{+}dt)\,\cdot(s)\rangle_c$", fontsize=10)
    fig.suptitle(r"$F^2$--$\sigma^2$ cross (CORRECTED, triple-sub, linear)  %s L1 %d cfg  op$_F$=%d"
                 % (tag, ncfg, OP_F), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_sigma2_cross_v2_lin_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)

    # log |C|
    figL, axsL = plt.subplots(1, 3, figsize=(14, 4.6))
    for ax, cc, lab, col, mk in zip(axsL, panels, LABELS,
                                    ["tab:red", "tab:green", "black"], ["o", "^", "*"]):
        ax.errorbar(dts, np.abs(cc[0][dts]), yerr=cc[1][dts], color=col, marker=mk, ms=5, lw=1, capsize=2)
        ax.set_yscale("log")
        ax.set_title(lab, fontsize=11)
        ax.set_xlabel(r"$dt$", fontsize=9)
        ax.grid(alpha=0.3, which="both")
    axsL[0].set_ylabel(r"$|\langle F^2(s{+}dt)\,\cdot(s)\rangle_c|$", fontsize=10)
    figL.suptitle(r"$F^2$--$\sigma^2$ cross (CORRECTED, triple-sub, LOG $|C|$)  %s L1 %d cfg  op$_F$=%d"
                  % (tag, ncfg, OP_F), fontsize=12)
    figL.tight_layout(rect=[0, 0, 1, 0.94])
    outL = "figs/f2_sigma2_cross_v2_log_%s_claude.png" % tag
    figL.savefig(outL, dpi=130)
    plt.close(figL)
    print("\n# -> %s" % out)
    print("# -> %s" % outL)


if __name__ == "__main__":
    main()
