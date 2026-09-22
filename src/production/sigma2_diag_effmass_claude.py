#!/usr/bin/env python3
# sigma2_diag_effmass_claude.py -- plain log-effmass of the DIAGONAL flavor x geometry sigma^2 correlators (no GEVP),
#   from the FULL 9-op nsrc2 cache.  Shows what each single operator {PP,FF} x {s2,O2m,O1m} overlaps before mixing.
#   Reads the _mc1 cache when MODE_CONTACT=1 (the contact fix); jackknife errors over BINSIZE bins.
#   Run: MODE_CONTACT=1 ENS=<L2> LREF=2 NVDIR=distill_Nv24_v2 BINSIZE=10 python3 sigma2_diag_effmass_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import glob
import re
import numpy as np
import distill_contract_claude as dc

SPLIT = int(os.environ.get("SPLIT", "1"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))
# 6 parity-even ops: flavor*3+geom, flavor 0=PP 1=FF ; geom 0=s2 1=O2m 2=O1m
OPS = [int(x) for x in os.environ.get("OPS", "0,1,2,3,4,5").split(",")]
OPLAB = {0: "PP-s2", 1: "PP-O2m", 2: "PP-O1m", 3: "FF-s2", 4: "FF-O2m", 5: "FF-O1m"}
# color-blind: distinct marker AND color per geometry; PP solid, FF open
COL = {0: "tab:blue", 1: "tab:red", 2: "tab:green", 3: "tab:blue", 4: "tab:red", 5: "tab:green"}
MRK = {0: "o", 1: "s", 2: "^", 3: "o", 4: "s", 5: "^"}
FILL = {0: "full", 1: "full", 2: "full", 3: "none", 4: "none", 5: "none"}


def jk_effmass(diag_cfg):
    # diag_cfg: (ncfg, T).  Bin, jackknife the log-effmass log(C(t)/C(t+1)).
    ncfg = diag_cfg.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([diag_cfg[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    C = blk.mean(0)
    with np.errstate(all="ignore"):
        em_c = np.log(C[:-1] / C[1:])
        ems = np.array([np.log(np.delete(blk, i, 0).mean(0)[:-1] / np.delete(blk, i, 0).mean(0)[1:]) for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    return em_c, em_e


def main():
    tag = dc.ENS.split("nu0")[0]
    CACHEDIR = "sigma2_flavor_cache_claude"
    MCTAG = "_mc1" if int(os.environ.get("MODE_CONTACT", "0")) else ""
    hits = glob.glob("%s/sigma2_flavorgeom_FULL_%s_L%d_*cfg_nsrc2_d%d%s_claude.npy" % (CACHEDIR, tag.replace(".", "p"), dc.L, SPLIT, MCTAG))
    if not hits:
        print("# no cache (MODE_CONTACT=%s MCTAG=%s); build sigma2_flavorgeom_full_v2_nsrc2 first" % (os.environ.get("MODE_CONTACT", "0"), MCTAG))
        return
    cache = max(hits, key=lambda p: int(re.search(r"_(\d+)cfg", p).group(1)))
    allC = np.load(cache)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    ncfg = allC.shape[0]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# loaded <- %s  shape %s   m_PS=%.4f 2m_PS=%.4f" % (cache, allC.shape, mps, m2ps))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.8)
    ax.text(TMAXPLOT * 0.72, m2ps + 0.012, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="dimgray", ls=":", lw=1.1, alpha=0.8)
    ax.text(TMAXPLOT * 0.72, mps + 0.012, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="dimgray")

    print("\n# diagonal log-effmass a_t m_eff(t):")
    for op in OPS:
        em_c, em_e = jk_effmass(allC[:, op, op, :].real)
        ts = np.arange(len(em_c))
        g = np.isfinite(em_c) & np.isfinite(em_e) & (em_e < 0.25)
        mf = "none" if FILL[op] == "none" else COL[op]
        ax.errorbar(ts[g] + 0.03 * op, em_c[g], yerr=em_e[g], color=COL[op], marker=MRK[op], ms=5, lw=1.0,
                    capsize=2, mfc=mf, label=OPLAB[op])
        plateau = "  ".join("%.3f" % em_c[t] for t in range(4, min(11, len(em_c))))
        print("#  %-8s dt4..10: %s" % (OPLAB[op], plateau))

    ax.set_ylim(0.2, 1.2)
    ax.set_xlim(1, TMAXPLOT)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"Diagonal $\sigma^2$ correlator effmass (NO GEVP)  %s L%d  MODE_CONTACT=%s  %dcfg"
                 % (tag, dc.L, os.environ.get("MODE_CONTACT", "0"), ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right", ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_diag_effmass_%s_L%d_mc%s_claude.png" % (tag, dc.L, os.environ.get("MODE_CONTACT", "0"))
    fig.savefig(out, dpi=140)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
