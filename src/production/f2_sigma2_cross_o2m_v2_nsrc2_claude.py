#!/usr/bin/env python3
# f2_sigma2_cross_o2m_v2_nsrc2_claude.py   [nsrc=2 DATA variant of f2_sigma2_cross_o2m_v2_claude.py]
#   Same corrected analysis (folded, vac-sub, triple subtraction) but on the nsrc=2 perambulators
#   (NVDIR=distill_Nv24_v2, source windows tsrc_list=[0,64]).  BOTH windows are used with the averaging
#   done STRICTLY POST-CONTRACTION: the equal-time antipodal loop L_2m(s) is built and the F^2 cross
#   contracted SEPARATELY within each window w (each already source-averages over s in that window), and
#   ONLY the finished per-window correlators C_w(dt) are averaged at fixed separation
#       C(dt) = (1/n_win) sum_w C_w(dt).
#   tau is NEVER combined across windows (make_config_win builds AblkS from that window's tau_w alone).
#   For nsrc=1 files (v1) load_peram_windows returns 1 window -> IDENTICAL to f2_sigma2_cross_o2m_v2 (the
#   machine-precision validation:  run this with NVDIR=distill_Nv24 and it must reproduce v1 exactly).
#
#   Env: WIN_ONLY=-1 (all windows; >=0 restricts to that single window, for a v2-window0-vs-v1 check).
#   Run: NVDIR=distill_Nv24_v2 python3 f2_sigma2_cross_o2m_v2_nsrc2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24_v2")
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
WIN_ONLY = int(os.environ.get("WIN_ONLY", "-1"))     # -1 = all windows ; >=0 = only that window


def gluedir():
    return "data_" + dc.ENS


def glue_has_O(k):
    # the per-timeslice F^2 density O (n_shapes, Nt) exists only on the config indices where the FULL glue
    # measurement was run; some v2-peram indices have a lighter glue file (F only, no O) -> must be skipped.
    p = "%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k)
    if not os.path.exists(p):
        return False
    try:
        with h5py.File(p, "r") as g:
            return "O" in g
    except Exception:
        return False


def load_of(k):
    with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (gluedir(), k), "r") as g:
        return np.array(g["O"])[OP_F].astype(float)     # (Nt,)


def o2m_loop_win(k, w, dual, Pmap):
    # per-window connected O_2m self-loop L_2m(s) (equal-time antipodal), built from window w's tau alone
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, w)
    idx = np.arange(nsite)
    L = np.full(twin, np.nan)
    for s in range(twin):
        Pf = AblkS(s, s)                                 # (nsite,NS,nsite,NS) equal-time
        B1 = Pf[idx, :, Pmap, :]                          # (nsite,NS,NS): P(x ; P(x))
        B2 = Pf[Pmap, :, idx, :]                          # (nsite,NS,NS): P(P(x) ; x)
        L[s] = -(dual * np.einsum('xab,xba->x', B1, B2)).sum().real
    return L, tsrc0, twin


def nwindows(k):
    V, windows = dc.load_peram_windows(k)
    return len(windows)


def cross_one_config(k, dual, Pmap):
    # POST-CONTRACTION window average: contract the F^2 x O_2m cross fully within each window, then mean.
    of = load_of(k)
    nw = nwindows(k)
    ws = [WIN_ONLY] if WIN_ONLY >= 0 else list(range(nw))
    acc = np.zeros(DTCALC)
    for w in ws:
        L, tsrc0, twin = o2m_loop_win(k, w, dual, Pmap)
        svalid = np.arange(twin)
        Cd = np.zeros(DTCALC)
        for dt in range(DTCALC):
            idxf = (tsrc0 + svalid + dt) % Nt
            idxb = (tsrc0 + svalid - dt) % Nt
            ofd = 0.5 * (of[idxf] + of[idxb])            # F^2 folded around THIS window's tsrc0
            Cd[dt] = np.mean(ofd * L[svalid])
        acc += Cd
    return acc / len(ws)


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh().astype(float)
    Pmap = G.antipodal_map()
    ks = [k for k in dc.KS if glue_has_O(k)]
    winlab = "all" if WIN_ONLY < 0 else "w%d" % WIN_ONLY
    print("# ENS=%s NVDIR=%s matched cfg=%d/%d OP_F=%d WIN=%s  (F^2 x antipodal O_2m loop, nsrc2)"
          % (tag, os.environ["NVDIR"], len(ks), len(dc.KS), OP_F, winlab))

    allC = []
    for k in ks:
        allC.append(cross_one_config(k, dual, Pmap))
    allC = np.array(allC)                                 # (ncfg, DTCALC)
    ncfg = allC.shape[0]

    ts = allC - allC[:, TSUM_LO:].mean(1, keepdims=True)
    n = ts.shape[0]
    samp = np.array([np.delete(ts, i, 0).mean(0) for i in range(n)])
    samp = samp - samp[:, PLAT_LO:DTCALC].mean(1, keepdims=True)
    cm = samp.mean(0)
    ee = np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    print("\n#  dt |   C(err)             S/N")
    for dt in range(0, DTMAX):
        sn = cm[dt] / ee[dt] if ee[dt] > 0 else 0.0
        print("#  %2d | % .3e(%.1e)  %6.2f" % (dt, cm[dt], ee[dt], sn))
    print("# max|S/N| = %.2f" % np.max(np.abs(cm[:DTMAX] / ee[:DTMAX])))

    at = float(os.environ.get("AT", "0.2"))
    with np.errstate(all="ignore"):
        em_samp = np.log(samp[:, :-1] / samp[:, 1:])
    em_c = em_samp.mean(0)
    em_e = np.sqrt((n - 1) * np.mean((em_samp - em_c) ** 2, 0))
    print("\n#  dt | a_t m_eff(err)   ->  m_phys=a_t m/a_t (err)   [F2 glue: a_t m=0.616 lat / 3.08 phys]")
    for dt in range(0, 12):
        if np.isfinite(em_c[dt]) and np.isfinite(em_e[dt]):
            print("#  %2d | %7.4f(%.4f)      ->  %6.3f(%.3f)"
                  % (dt, em_c[dt], em_e[dt], em_c[dt] / at, em_e[dt] / at))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(0, DTMAX)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(dts, np.abs(cm[dts]), yerr=ee[dts], color="tab:blue", marker="s", ms=5, lw=1, capsize=2)
    ax.set_yscale("log")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$|\langle F^2(s{+}dt)\,O_{2m}(s)\rangle_c|$")
    ax.grid(alpha=0.3)
    ax.set_title(r"$F^2$ x antipodal $O_{2m}$, folded, vac-sub, nsrc2 (%s)  %s L1 %d cfg"
                 % (winlab, tag, ncfg), fontsize=11)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_sigma2_cross_o2m_v2_nsrc2_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
