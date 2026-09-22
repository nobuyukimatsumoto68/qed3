#!/usr/bin/env python3
# sigma2_Peven_partialhankel_scan_claude.py -- 5x5 scan of PER-OPERATOR Hankel offset-SETS for the two ops
#   that best split the two-meson: s2 (=sigma^2_00) and O1m (coincident single-sum).  Each op independently
#   takes one of 5 offset-sets {0-2, 0-3, 0-4, 0-2-4, 0-3-6}; the op is replicated at each NONZERO offset in
#   its set (off=0 copy is already the base op), both flavors PP,FF.  25 combinations -> 25 GEVP effmass
#   panels in ONE png (rows = s2 set, cols = O1m set).  O2m and the base copies stay plain (off 0).
#   Big(t)[a,b] = C_base[op_a,op_b](t + off_a + off_b) ; then plain rebased GEVP (reb NKEEP @ REBT, T0).
#   Run: ENS=.. NVDIR=.. LREF=2 REBT=4 NKEEP=3 T0=3 BINSIZE=10 TMAXPLOT=14 python3 sigma2_Peven_partialhankel_scan_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24_v2")
import sys
sys.path.insert(0, ".")
import glob
import re
import numpy as np
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs

REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "3"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
SPLIT = int(os.environ.get("SPLIT", "1"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "14"))

# cache op index = flavor*3 + geom ; geom 0=s2 1=O2m 2=O1m ; flavor 0=PP 1=FF.
IDX_S2 = [0, 3]                                             # PP-s2, FF-s2
IDX_O1M = [2, 5]                                            # PP-O1m, FF-O1m
BASE = [0, 1, 2, 3, 4, 5]                                   # P+ {PP,FF} x {s2,O2m,O1m}, all off 0
OFFSET_SETS = {"0-2": [0, 2], "0-3": [0, 3], "0-4": [0, 4], "0-2-4": [0, 2, 4], "0-3-6": [0, 3, 6]}
ALL_KEYS = ["0-2", "0-3", "0-4", "0-2-4", "0-3-6"]
S2_KEYS = os.environ.get("S2_KEYS", ",".join(ALL_KEYS)).split(",")     # rows ; comma-list subset of ALL_KEYS
O1M_KEYS = os.environ.get("O1M_KEYS", ",".join(ALL_KEYS)).split(",")   # cols


def build_basis(set_s2, set_o1m):
    baseidx = list(BASE)
    off = [0] * 6
    for o in set_s2:
        if o == 0:
            continue
        baseidx += IDX_S2                                   # s2 replica at offset o (both flavors)
        off += [o, o]
    for o in set_o1m:
        if o == 0:
            continue
        baseidx += IDX_O1M                                  # O1m replica at offset o
        off += [o, o]
    return baseidx, off


def build_big(Cmat, baseidx, off):
    n = len(baseidx)
    omax = max(off)
    DT = Cmat.shape[2]
    tmax = DT - 2 * omax
    Big = np.full((tmax, n, n), np.nan)
    for t in range(tmax):
        for a in range(n):
            for b in range(n):
                Big[t, a, b] = Cmat[baseidx[a], baseidx[b], t + off[a] + off[b]]
    return Big


def gevp(Cmat, baseidx, off, Vfix):
    Big = build_big(Cmat, baseidx, off)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def scan_one(blk, baseidx, off):
    nb = blk.shape[0]
    em_c, Vfix = gevp(blk.mean(0), baseidx, off, None)
    ems = np.array([gevp(np.delete(blk, i, 0).mean(0), baseidx, off, Vfix)[0] for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    return em_c, em_e


def main():
    tag = dc.ENS.split("nu0")[0]
    CACHEDIR = "sigma2_flavor_cache_claude"
    MCTAG = "_mc1" if int(os.environ.get("MODE_CONTACT", "0")) else ""   # mode-space contact fix -> _mc1 cache
    hits = glob.glob("%s/sigma2_flavorgeom_FULL_%s_L%d_*cfg_nsrc2_d%d%s_claude.npy" % (CACHEDIR, tag.replace(".", "p"), dc.L, SPLIT, MCTAG))
    if not hits:
        hits = glob.glob("%s/sigma2_flavorgeom_FULL_%s_*cfg_d%d%s_claude.npy" % (CACHEDIR, tag.replace(".", "p"), SPLIT, MCTAG))
    cache = max(hits, key=lambda p: int(re.search(r"_(\d+)cfg", p).group(1)))
    allC = np.load(cache)
    print("# loaded <- %s  shape %s" % (cache, allC.shape))
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m2ps = float(os.environ.get("M2PS", "0.7054"))         # 2 m_PS ; L2 Nf2 g1 = 0.7054(28) ; L1 = 0.644
    cols = ["tab:green", "tab:red", "tab:blue"]
    mkr = ["o", "s", "^"]
    nr = len(S2_KEYS)
    nc = len(O1M_KEYS)
    fig, axes = plt.subplots(nr, nc, figsize=(4.2 * nc, 3.6 * nr), sharex=True, sharey=True, squeeze=False)
    for i, ks2 in enumerate(S2_KEYS):
        for j, ko1m in enumerate(O1M_KEYS):
            ax = axes[i, j]
            baseidx, off = build_basis(OFFSET_SETS[ks2], OFFSET_SETS[ko1m])
            em_c, em_e = scan_one(blk, baseidx, off)
            tmax = em_c.shape[0]
            ts = np.arange(tmax)
            ax.axhline(m2ps, color="gray", ls="--", lw=1, alpha=0.6)
            for n in range(NKEEP):
                g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
                ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n], marker=mkr[n], ms=4, lw=1.0,
                            capsize=2.0, label="state %d" % n if (i == 0 and j == 0) else None)
            ax.axvline(REBT, color="k", ls=":", lw=0.8, alpha=0.3)
            ax.set_ylim(0.25, 1.05)
            ax.set_xlim(T0, min(tmax, TMAXPLOT))
            ax.set_title("s2:%s  O1m:%s  (%dop)" % (ks2, ko1m, len(baseidx)), fontsize=9)
            ax.grid(alpha=0.3)
            if j == 0:
                ax.set_ylabel(r"$a_t m_\mathrm{eff}$", fontsize=9)
            if i == nr - 1:
                ax.set_xlabel(r"$t$", fontsize=9)
            print("# s2:%-5s O1m:%-5s  m0..%d @t=8: %s" % (ks2, ko1m, NKEEP - 1,
                  "  ".join("%.4f(%.4f)" % (em_c[8, n], em_e[8, n]) if 8 < em_c.shape[0] else "nan" for n in range(NKEEP))))
    axes[0, 0].legend(fontsize=8, loc="upper right")
    fig.suptitle(r"P+ 6op + per-op Hankel offset-set scan (s2 rows, O1m cols)  reb%d@%d T0=%d  %s L%d %dcfg  ($2m_{PS}=%.3f$)"
                 % (NKEEP, REBT, T0, tag, dc.L, ncfg, m2ps), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    os.makedirs("figs", exist_ok=True)
    subtag = "" if (len(S2_KEYS) == 5 and len(O1M_KEYS) == 5) else "_s2-%s_o1m-%s" % ("_".join(k.replace("-", "") for k in S2_KEYS), "_".join(k.replace("-", "") for k in O1M_KEYS))
    out = "figs/sigma2_Peven_partialhankel_scan%s_%s_L%d_reb%d_T0%d_claude.png" % (subtag, tag, dc.L, NKEEP, T0)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
