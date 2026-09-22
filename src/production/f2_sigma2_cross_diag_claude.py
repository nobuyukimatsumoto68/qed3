#!/usr/bin/env python3
# f2_sigma2_cross_diag_claude.py
#   Diagram-by-diagram F^2 -- sigma^2 CROSS correlator  <F^2(s+dt) sigma^2(s)>, distillation, PS channel.
#   F^2 is purely gluonic, so a single sigma^2(s) contracts into TWO fermionic topologies (linked to
#   F^2 only through the gauge field):
#       disc  :  D_S(s)^2      (two separate sigma tadpole loops ; contact-subtracted -> <D_S>=0)
#       ext   :  D'_S(s)       (one connected double-loop)
#   with  D_S(s)  = Tr[Phi(s) tt(s,s)] ,  D'_S(s) = Tr[Phi(s) tt(s,s) Phi(s) tt(s,s)] ,
#         tt(s,s) = tau(s,s) - CONTACT*I   (GW contact 1/2 ; CONTACT=0.5).
#   sigma^2 total (PS, one_point convention) = 2 ( D_S^2 + D'_S ).
#
#   O_F(t) = glue F^2 shape operator OP_F (l=0, p2, op 0 = basic F^2), full Nt.  Translation-average the
#   source s over the perambulator window [tsrc0, tsrc0+twin); F^2 at (tsrc0+s+dt) % Nt.
#   t-sum -> plateau subtracted per jackknife ensemble (removes the <O_F><X> vacuum constant); LINEAR.
#
#   Run: ENS=Nf2_gsq1.000000...L1_hb1.000000 NVDIR=distill_Nv24 CONTACT=0.5 python3 f2_sigma2_cross_diag_claude.py

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
DTMAX = int(os.environ.get("DTMAX", "24"))          # plotted range
DTCALC = int(os.environ.get("DTCALC", "48"))        # computed range (for the plateau tail)
PLAT_LO = int(os.environ.get("PLAT_LO", "32"))      # plateau (t-sum) window start


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
    # cross per dt: translation-average source a over window
    Cd = np.zeros((2, DTCALC))                             # [disc, ext]
    for dt in range(DTCALC):
        idx = (tsrc0 + np.arange(twin) + dt) % Nt
        ofd = of[idx]
        Cd[0, dt] = np.mean(ofd * (DS ** 2))              # disc  D_S^2
        Cd[1, dt] = np.mean(ofd * DpS)                    # ext   D'_S
    return Cd


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    print("# ENS=%s  matched cfg=%d/%d  OP_F=%d CONTACT=%.2f" % (tag, len(ks), len(dc.KS), OP_F, CONTACT))
    allC = np.array([per_config(k, w00) for k in ks])     # (ncfg, 2, DTCALC)
    ncfg = allC.shape[0]

    def jk_corr(C):
        # C: (ncfg, DTCALC).  plateau (t-sum) subtraction + config jackknife.
        plat = C[:, PLAT_LO:].mean(1, keepdims=True)
        Cc = C - plat
        n = C.shape[0]
        samp = np.array([np.delete(Cc, i, 0).mean(0) for i in range(n)])
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    disc = jk_corr(allC[:, 0, :])
    ext = jk_corr(allC[:, 1, :])
    tot = jk_corr(2.0 * (allC[:, 0, :] + allC[:, 1, :]))  # sigma^2 total = 2(D_S^2 + D'_S)

    print("\n#  dt |   disc(D_S^2)        ext(D'_S)         TOTAL 2(disc+ext)")
    for dt in range(0, DTMAX):
        print("#  %2d | %11.3e(%.1e) %11.3e(%.1e) %11.3e(%.1e)"
              % (dt, disc[0][dt], disc[1][dt], ext[0][dt], ext[1][dt], tot[0][dt], tot[1][dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(0, DTMAX)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6))
    panels = [(disc, "disc  $D_S^2$", "tab:red", "o"),
              (ext, "ext  $D'_S$", "tab:green", "^"),
              (tot, "TOTAL  $2(D_S^2+D'_S)$", "black", "*")]
    for ax, (cc, lab, col, mk) in zip(axs, panels):
        ax.errorbar(dts, cc[0][dts], yerr=cc[1][dts], color=col, marker=mk, ms=5, lw=1, capsize=2)
        ax.axhline(0.0, color="gray", lw=0.8, alpha=0.6)
        ax.set_title(lab, fontsize=11)
        ax.set_xlabel(r"$dt$", fontsize=9)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[0].set_ylabel(r"$\langle F^2(s{+}dt)\,\cdot(s)\rangle_c$", fontsize=10)
    fig.suptitle(r"$F^2$--$\sigma^2$ cross, per diagram (PS, linear, t-sum+plateau sub)  %s L1 %d cfg  op$_F$=%d"
                 % (tag, ncfg, OP_F), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_sigma2_cross_diag_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
