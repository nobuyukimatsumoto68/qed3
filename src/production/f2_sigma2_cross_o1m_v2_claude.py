#!/usr/bin/env python3
# f2_sigma2_cross_o1m_v2_claude.py
#   Cross  <F^2(s+dt) O_1m(s)>_c  with the EXISTING coincident time-split one-meson interpolator
#       O_1m(s) = sum_x A_x sigma(x,s) sigma(x,s+delta)   (fs_channels_v2 / fs_gevp_point op 2).
#   F^2 is gluonic -> O_1m's four fermions self-contract into a single connected loop (the disc tadpole
#   product dies under improvement).  The connected self-loop (one number per source s):
#       L_delta(s) = - sum_x A_x Tr_spin[ P(x,s ; x,s+delta) P(x,s+delta ; x,s) ]
#   with the position propagator block P = AblkS(ta,tb) = U_ta tau(ta,tb) U_tb^dag (improved on diagonal;
#   here ta != tb so no contact).  This is the SAME O_1m loop as the point-operator machinery, opened up
#   in time by delta so the fermion loop propagates -> filters toward the lightest 0++ meson content
#   BEFORE asking whether F^2 couples to it.  (F^2 cannot reach the 4-fermion two-meson STATE by particle
#   number; this probes the gluonic admixture of the light single-loop 0++ density.)
#
#   Clean t-fold (C(dt)=C(-dt), F^2 global).  Per-config t-sum + jackknife plateau (vacuum) subtraction.
#   Run: NVDIR=distill_Nv24 python3 f2_sigma2_cross_o1m_v2_claude.py     (SPLITS=1,2,4 default)

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
import fs_gevp_point_claude as G

Nt = 128
OP_F = int(os.environ.get("OP_F", "0"))
DTCALC = int(os.environ.get("DTCALC", "56"))
DTMAX = int(os.environ.get("DTMAX", "40"))
TSUM_LO = int(os.environ.get("TSUM_LO", "1"))
PLAT_LO = int(os.environ.get("PLAT_LO", str(DTCALC - 8)))
SPLITS = [int(x) for x in os.environ.get("SPLITS", "1,2,4").split(",")]


def gluedir():
    return "data_" + dc.ENS


def load_of(k):
    with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k), "r") as g:
        return np.array(g["O"])[OP_F].astype(float)     # (Nt,)


def o1m_loops(k, dual):
    # per-source connected O_1m self-loop L_delta(s) for each delta in SPLITS ; also tsrc0, twin
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(k)
    idx = np.arange(nsite)
    V, _, _, tsrc0, _ = dc.load_peram(k)                # tsrc0 for the F^2 index
    loops = {}
    for d in SPLITS:
        L = np.full(twin, np.nan)
        for s in range(twin - d):
            Pf = AblkS(s, s + d)                         # (nsite,NS,nsite,NS)
            Pb = AblkS(s + d, s)
            Pf_d = Pf[idx, :, idx, :]                    # (nsite,NS,NS) site-coincident block
            Pb_d = Pb[idx, :, idx, :]
            L[s] = -(dual * np.einsum('xab,xba->x', Pf_d, Pb_d)).sum().real
        loops[d] = L
    return loops, tsrc0, twin


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh().astype(float)
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k))]
    print("# ENS=%s matched cfg=%d/%d OP_F=%d SPLITS=%s  (F^2 x coincident-split O_1m loop)"
          % (tag, len(ks), len(dc.KS), OP_F, SPLITS))

    allC = {d: [] for d in SPLITS}                        # per delta: list over cfg of (DTCALC,)
    for k in ks:
        loops, tsrc0, twin = o1m_loops(k, dual)
        of = load_of(k)
        for d in SPLITS:
            L = loops[d]
            svalid = np.array([s for s in range(twin - d) if np.isfinite(L[s])])
            Cd = np.zeros(DTCALC)
            for dt in range(DTCALC):
                idxf = (tsrc0 + svalid + dt) % Nt
                idxb = (tsrc0 + svalid - dt) % Nt
                ofd = 0.5 * (of[idxf] + of[idxb])         # folded F^2 around the first sigma (s)
                Cd[dt] = np.mean(ofd * L[svalid])
            allC[d].append(Cd)
    for d in SPLITS:
        allC[d] = np.array(allC[d])                       # (ncfg, DTCALC)
    ncfg = allC[SPLITS[0]].shape[0]

    def sub_jk(C):
        # per-config t-sum sub -> jackknife + plateau (vacuum tail) sub
        ts = C - C[:, TSUM_LO:].mean(1, keepdims=True)
        n = ts.shape[0]
        samp = np.array([np.delete(ts, i, 0).mean(0) for i in range(n)])
        samp = samp - samp[:, PLAT_LO:DTCALC].mean(1, keepdims=True)
        return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    res = {d: sub_jk(allC[d]) for d in SPLITS}

    print("\n# ext(O_1m) cross  S/N by dt (per delta):")
    hdr = "#  dt |" + "".join("  d=%d  C(err)          S/N" % d for d in SPLITS)
    print(hdr)
    for dt in range(0, DTMAX):
        row = "#  %2d |" % dt
        for d in SPLITS:
            cm, ee = res[d]
            sn = cm[dt] / ee[dt] if ee[dt] > 0 else 0.0
            row += "  % .3e(%.1e) %6.2f" % (cm[dt], ee[dt], sn)
        print(row)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(0, DTMAX)
    cols = {1: "tab:red", 2: "tab:blue", 4: "tab:green", 8: "tab:purple"}
    mks = {1: "o", 2: "s", 4: "^", 8: "D"}
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    for d in SPLITS:
        cm, ee = res[d]
        c = cols.get(d, "black")
        m = mks.get(d, "o")
        axs[0].errorbar(dts, cm[dts], yerr=ee[dts], color=c, marker=m, ms=5, lw=1, capsize=2, label=r"$\delta=%d$" % d)
        axs[1].errorbar(dts, np.abs(cm[dts]), yerr=ee[dts], color=c, marker=m, ms=5, lw=1, capsize=2, label=r"$\delta=%d$" % d)
    axs[0].axhline(0.0, color="gray", lw=0.8, alpha=0.6)
    axs[0].set_title("linear", fontsize=11)
    axs[1].set_yscale("log")
    axs[1].set_title(r"$|C|$ log", fontsize=11)
    for ax in axs:
        ax.set_xlabel(r"$dt$")
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel(r"$\langle F^2(s{+}dt)\,O_{1m}(s)\rangle_c$")
    fig.suptitle(r"$F^2$ x coincident-split $O_{1m}$ (one-meson interpolator), folded, vac-sub  %s L1 %d cfg"
                 % (tag, ncfg), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_sigma2_cross_o1m_v2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
