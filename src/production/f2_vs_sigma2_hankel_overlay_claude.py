#!/usr/bin/env python3
# f2_vs_sigma2_hankel_overlay_claude.py
#   Overlay the GLUONIC F^2 glueball on the trusted {sigma^2_00, O_2m, O_1m} Hankel+rebase GEVP curves.
#   Same Hankel machinery as fs_channels_v2_hankel_claude.py (OFFSETS=0,2,4 REBT=4 NKEEP=2 T0=3).
#   Units: LATTICE a_t m (the fermionic native units).  F^2 glueball = 2.88(12) physical -> x a_t = 0.576(24)
#   lattice.  Shows where the heavy gluonic 0++ sits vs the sigma^2 one-meson-rich ground (~0.46) and the
#   two-meson (~0.62 ~ 2 m_PS).
#   Run: python3 f2_vs_sigma2_hankel_overlay_claude.py [glue_dat]

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import glob
import numpy as np
import h5py
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs
import fs_gevp_point_claude as G

AT = float(os.environ.get("AT", "0.2"))
T0 = int(os.environ.get("T0", "3"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "2"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NOP = int(os.environ.get("NOP", "3"))
SPLIT = int(os.environ.get("SPLIT", "1"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,2,4").split(",")]
M2PS = 0.644                                       # lattice 2 m_PS reference
F2_PHYS = 2.8816
F2_PHYS_E = 0.1178
GLUE_DAT = sys.argv[1] if len(sys.argv) > 1 else \
    "/tmp/claude-1000/-mnt-barracuda22-qed3/cdb76ded-4ed0-4f78-8eda-57542aa97f37/scratchpad/f2_gevp_repro.dat"


def hankel_reb(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def cross_effmass_lattice(tag):
    # F^2 x O_2m cross correlator effmass (LATTICE a_t m), folded + per-config t-sum + plateau vacuum-sub
    Nt = 128
    DTCALC = 48
    ddir = "data_" + dc.ENS
    dual = dc.dual_areas_from_mesh().astype(float)
    Pmap = G.antipodal_map()
    ks = [k for k in dc.KS if os.path.exists("%s/glue_f2_v2_shapes.%d.h5" % (ddir, k))]
    allC = np.zeros((len(ks), DTCALC))
    for ic, k in enumerate(ks):
        with h5py.File("%s/glue_f2_v2_shapes.%d.h5" % (ddir, k), "r") as g:
            of = np.array(g["O"])[0].astype(float)
        AblkS, _, twin, nsite, U, tau = G.make_config(k)
        idx = np.arange(nsite)
        V, _, _, tsrc0, _ = dc.load_peram(k)
        L = np.empty(twin)
        for s in range(twin):
            Pf = AblkS(s, s)
            L[s] = -(dual * np.einsum('xab,xba->x', Pf[idx, :, Pmap, :], Pf[Pmap, :, idx, :])).sum().real
        sv = np.arange(twin)
        for dt in range(DTCALC):
            idxf = (tsrc0 + sv + dt) % Nt
            idxb = (tsrc0 + sv - dt) % Nt
            allC[ic, dt] = np.mean(0.5 * (of[idxf] + of[idxb]) * L)
    ts = allC - allC[:, 1:].mean(1, keepdims=True)               # per-config t-sum sub
    nb = len(ks) // BINSIZE
    blk = np.array([ts[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    samp = np.array([np.delete(blk, i, 0).mean(0) for i in range(nb)])
    samp = samp - samp[:, DTCALC - 8:].mean(1, keepdims=True)    # plateau sub
    with np.errstate(all="ignore"):
        em = np.log(samp[:, :-1] / samp[:, 1:])                  # both negative & decaying -> real
    return em.mean(0), np.sqrt((nb - 1) * np.mean((em - em.mean(0)) ** 2, 0))


def main():
    tag = dc.ENS.split("nu0")[0]
    cands = glob.glob("fs_channels_v2_cache_claude/fs_channels_v2_%s_*cfg_d%d_claude.npy"
                      % (tag.replace(".", "p"), SPLIT))
    cache = max(cands, key=os.path.getsize)
    allC = np.load(cache)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))[:, :NOP, :NOP, :]
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = hankel_reb(blk.mean(0), None)
    ems = np.array([hankel_reb(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]
    print("# {sigma^2_00,O_2m,O_1m} Hankel+rebase (lattice a_t m):  F^2 glueball = %.3f(%.3f) lattice"
          % (F2_PHYS * AT, F2_PHYS_E * AT))
    for t in range(1, min(tmax, 20)):
        print("#  t=%2d | " % t + "  ".join("m%d=%.4f(%.4f)" % (n, em_c[t, n], em_err[t, n]) for n in range(NKEEP)))

    # F^2 glueball effmass (physical from glue GEVP .dat state 0) -> lattice x a_t
    g = np.loadtxt(GLUE_DAT)
    gt = np.rint(g[:, 0] / AT).astype(int)          # timeslice
    gF = g[:, 3] * AT                               # physical -> lattice
    gFe = g[:, 4] * AT

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    # 2 m_PS and F^2 references
    ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(16.5, M2PS + 0.006, r"$2m_{PS}\approx0.644$", fontsize=8.5, color="gray")
    ax.axhspan(F2_PHYS * AT - F2_PHYS_E * AT, F2_PHYS * AT + F2_PHYS_E * AT, color="black", alpha=0.10)
    ax.axhline(F2_PHYS * AT, color="black", ls="-.", lw=1.1, alpha=0.7)
    ax.text(16.5, F2_PHYS * AT + 0.006, r"$F^2$ fit $0.576$", fontsize=8.5, color="black")
    cols = ["tab:green", "tab:red"]
    mkr = ["o", "s"]
    labs = [r"$\sigma^2$ state 0 (1-meson-rich)", r"$\sigma^2$ state 1 (two-meson)"]
    for n in range(NKEEP):
        gmask = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.25)
        ax.errorbar(ts[gmask], em_c[gmask, n], yerr=em_err[gmask, n], color=cols[n], marker=mkr[n],
                    ms=5, lw=1.1, capsize=2.5, label=labs[n])
    gmask = (gFe < 0.20) & (gt >= 1) & (gt <= 8)
    ax.errorbar(gt[gmask], gF[gmask], yerr=gFe[gmask], color="black", marker="D", ms=6, lw=1.4,
                capsize=3, label=r"$F^2$ glueball GEVP")
    # F^2 x O_2m cross-channel effmass (lattice)
    xm, xe = cross_effmass_lattice(tag)
    xt = np.arange(len(xm))
    cm = np.isfinite(xm) & np.isfinite(xe) & (xe < 0.35) & (xt <= 4)
    ax.errorbar(xt[cm], xm[cm], yerr=xe[cm], color="tab:purple", marker="v", ms=8, lw=1.4, capsize=3.5,
                label=r"$F^2\!\times\!O_{2m}$ cross")
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.2, 0.95)
    ax.set_xlim(-0.6, 20)
    ax.set_xlabel(r"$t$ (timeslice)")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$ (lattice)")
    ax.set_title(r"$F^2$ glueball vs $\{\sigma^2_{00},O_{2m},O_{1m}\}$ Hankel+rebase  Dt=%s reb %d@%d T0=%d  %s L1 %d cfg"
                 % (OFFSETS, NKEEP, REBT, T0, tag, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_vs_sigma2_hankel_overlay_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
