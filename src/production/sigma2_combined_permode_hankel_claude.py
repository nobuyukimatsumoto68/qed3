#!/usr/bin/env python3
# sigma2_combined_permode_hankel_claude.py -- combined (2,2)-shell + sigma^2 GEVP with PER-OPERATOR block-Hankel.
#   Each operator carries its own offset list (hs.hankel_permode); the raw 6x6xDTMAX per-config correlator
#   tensor allC is cached to .npy so Hankel-knob iteration is instant (the expensive cross-blob recompute runs once).
#   Basis order 0..5 = {shell l1/2, shell l3/2 (2,2), shell l5/2, s2, O2m, O1m}.
#   Refs: block-Hankel / GPOF Aubin-Orginos arXiv:1010.0202.
#   Run: MODE_CONTACT=1 ENS=<L2> LREF=2 NVDIR=distill_Nv24_v2 WINDOWS=0-4,4-12,12-24 OPS2=0,1,2 \
#        OFFS_LIST="0-4;0-3;0-2-4;0-2;0-2;0-2" NKEEP=4 REBT=3 T0=2 BINSIZE=10 NCFG=200 \
#        python3 sigma2_combined_permode_hankel_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MODE_CONTACT", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import hankel_rebase_scan_claude as hs
import sigma2_combined_gevp_claude as CB

WINDOWS = CB.WINDOWS
OPS2 = CB.OPS2
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
REBT = int(os.environ.get("REBT", "3"))
T0 = int(os.environ.get("T0", "2"))
NKEEP = int(os.environ.get("NKEEP", "4"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))
# per-operator offsets: groups separated by ';', offsets within a group by '-'
OFFS_LIST = [[int(x) for x in grp.split("-")] for grp in os.environ.get("OFFS_LIST", "0-4;0-3;0-2-4;0;0-2;0-2").split(";")]
OPNAMES = ["shell l1/2", "shell l3/2 (2,2)", "shell l5/2", "s2", "O2m", "O1m"]
# KEEP_OPS: subset of the cached 6-op basis to use (indices into OPNAMES); OFFS_LIST is then per KEPT op.
KEEP_OPS = [int(x) for x in os.environ.get("KEEP_OPS", "0,1,2,3,4,5").split(",")]


def load_or_compute_allC(ks, dualf, wY, Pmap):
    tag = dc.ENS.split("nu0")[0]
    wtag = "_".join("%d-%d" % w for w in WINDOWS)
    otag = "".join(str(o) for o in OPS2)
    os.makedirs("cache_claude", exist_ok=True)
    cf = "cache_claude/combined_allC_%s_L%d_win%s_ops%s_mc%s_n%d.npy" % (
        tag, dc.L, wtag, otag, os.environ.get("MODE_CONTACT"), len(ks))
    if os.path.exists(cf):
        print("# loading cached allC <- %s" % cf)
        return np.load(cf), cf
    print("# computing allC (expensive cross blob) for %d cfg ..." % len(ks))
    allC = np.array([CB.one_config(k, dualf, wY, Pmap) for k in ks])
    np.save(cf, allC)
    print("# cached allC -> %s" % cf)
    return allC, cf


def gevp_permode(Cmat):
    Cts = np.transpose(Cmat, (2, 0, 1))                 # (DTMAX, ntot, ntot)
    Big = hs.hankel_permode(Cts, OFFS_LIST)
    Vp = hs.staged_project(Big, [(REBT, NKEEP)], T0)
    return hs.rebased_effmass_fixed(Big, Vp, T0)


def main():
    tag = dc.ENS.split("nu0")[0]
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# ENS=%s L=%d  PER-MODE HANKEL combined GEVP  windows=%s ops2=%s  %d cfg  reb%d@%d T0=%d NKEEP=%d"
          % (tag, dc.L, WINDOWS, OPS2, len(ks), REBT, NKEEP, T0, NKEEP))
    assert len(OFFS_LIST) == len(KEEP_OPS), "OFFS_LIST needs one group per kept op"
    for i, io in enumerate(KEEP_OPS):
        print("#   %-16s offsets=%s" % (OPNAMES[io], OFFS_LIST[i]))
    allC, cf = load_or_compute_allC(ks, dualf, wY, Pmap)
    allC = allC[:, KEEP_OPS][:, :, KEEP_OPS]
    ncfg = allC.shape[0]
    nb = max(ncfg // BINSIZE, 1)
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c = gevp_permode(blk.mean(0))
    if nb < 2:
        em_e = np.zeros_like(em_c)
    else:
        ems = np.array([gevp_permode(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
        em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax, nk = em_c.shape
    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(nk)))
    for t in range(T0, min(tmax, 16)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(nk))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.7)
    ax.text(TMAXPLOT * 0.7, m2ps + 0.01, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="firebrick", ls=":", lw=1.2, alpha=0.8)
    ax.text(TMAXPLOT * 0.7, mps + 0.01, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="firebrick")
    cols = ["firebrick", "tab:purple", "tab:blue", "tab:green", "tab:orange"]
    mkr = ["o", "D", "s", "^", "v"]
    for n in range(nk):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 5], marker=mkr[n % 5], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.set_ylim(0.2, 1.05)
    ax.set_xlim(T0, TMAXPLOT)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"GEVP ops=%s offs=%s  %s L%d %dcfg  reb%d@%d T0=%d"
                 % ([OPNAMES[i] for i in KEEP_OPS], OFFS_LIST, tag, dc.L, ncfg, REBT, NKEEP, T0), fontsize=8)
    ax.legend(fontsize=9, loc="upper right", ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    otag = "ops" + "".join(str(o) for o in KEEP_OPS) + "_off" + "_".join("".join(str(x) for x in g) for g in OFFS_LIST)
    out = "figs/sigma2_combined_permode_hankel_%s_L%d_%s_nk%d_claude.png" % (tag, dc.L, otag, NKEEP)
    fig.savefig(out, dpi=140)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
