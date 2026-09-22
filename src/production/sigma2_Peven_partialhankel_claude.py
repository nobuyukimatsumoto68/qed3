#!/usr/bin/env python3
# sigma2_Peven_partialhankel_claude.py -- P+ 6x6 GEVP with a PER-OPERATOR (partial) Hankel: only the ANTIPODAL
#   operator O_2m gets a time-shifted partner (offset HOFF), the others stay plain (offset 0).  Motivation:
#   resolve the TWO near-degenerate two-meson levels expected around 2 m_PS -- the antipodal op best overlaps
#   the two-meson, so a Hankel partner on IT adds the variational direction to split the pair, while the full
#   Hankel (all ops) destabilised the noisy L2 basis.  Enlarged basis = 6 base ops (off 0) + {PP-O2m, FF-O2m}
#   (off HOFF) = 8 ops.  Big(t)[a,b] = C_base[op_a,op_b](t + off_a + off_b).  Then plain rebased GEVP.
#   Run: ENS=.. NVDIR=.. LREF=2 HOFF=3 REBT=4 NKEEP=3 T0=3 BINSIZE=10 TMAXPLOT=18 python3 sigma2_Peven_partialhankel_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24_v2")
import sys
sys.path.insert(0, ".")
import glob
import re
import numpy as np
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs

HOFF = int(os.environ.get("HOFF", "3"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "3"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
SPLIT = int(os.environ.get("SPLIT", "1"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))

# enlarged basis: cache op index (flavor*3+geom) + per-op time offset.  geom: 0=s2,1=O2m,2=O1m ; flavor 0=PP 1=FF.
# HGEOM accepts a comma-list, e.g. "0,2" gives BOTH sigma_00 and O_1m a PP+FF off=HOFF Hankel partner.
HGEOMS = [int(x) for x in os.environ.get("HGEOM", "1").split(",")]   # geometries that get the Hankel partner
GLAB = ["s2", "O2m", "O1m"]
GLABSEL = "+".join(GLAB[hg] for hg in HGEOMS)
BASE = [0, 1, 2, 3, 4, 5]                                   # P+ {PP,FF} x {s2,O2m,O1m}
BASE_OFF = [0] * 6
HANKEL = []                                                 # per geom: PP-<geom>, FF-<geom> get the off=HOFF partner
for hg in HGEOMS:
    HANKEL.append(hg)
    HANKEL.append(hg + 3)
BASEIDX = BASE + HANKEL
LAB = ["PP-s2", "PP-O2m", "PP-O1m", "FF-s2", "FF-O2m", "FF-O1m"]
for hg in HGEOMS:
    LAB.append("PP-%s+%d" % (GLAB[hg], HOFF))
    LAB.append("FF-%s+%d" % (GLAB[hg], HOFF))


def build_big(Cmat, off):
    n = len(BASEIDX)
    omax = max(off)
    DT = Cmat.shape[2]
    tmax = DT - 2 * omax
    Big = np.full((tmax, n, n), np.nan)
    for t in range(tmax):
        for a in range(n):
            for b in range(n):
                Big[t, a, b] = Cmat[BASEIDX[a], BASEIDX[b], t + off[a] + off[b]]
    return Big


def gevp(Cmat, off, Vfix):
    Big = build_big(Cmat, off)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


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
    off = BASE_OFF + [HOFF] * len(HANKEL)
    print("# PARTIAL Hankel: %s get off=%d partner ; others off=0 ; %d ops ; reb%d@%d T0=%d"
          % (GLABSEL, HOFF, len(BASEIDX), NKEEP, REBT, T0))
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = gevp(blk.mean(0), off, None)
    ems = np.array([gevp(np.delete(blk, i, 0).mean(0), off, Vfix)[0] for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]
    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)))
    for t in range(T0, min(tmax, 20)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(NKEEP))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    m2ps = float(os.environ.get("M2PS", "0.7054"))     # 2 m_PS threshold; L2 Nf2 g1 = 0.7054(28) (m_PS=0.3527, scalar-5b) ; L1 = 0.644
    ax.axhline(m2ps, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(tmax * 0.6, m2ps + 0.008, r"$2m_{PS}=%.3f$ (L%d)" % (m2ps, dc.L), fontsize=8, color="gray")
    cols = ["tab:green", "tab:red", "tab:blue", "tab:orange", "tab:purple"]
    mkr = ["o", "s", "^", "D", "v"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 5], marker=mkr[n % 5], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.2, 1.1)
    ax.set_xlim(T0, min(tmax, TMAXPLOT))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"P+ 6op + %s Hankel(off=%d)  reb%d@%d T0=%d  %s L%d %dcfg"
                 % (GLABSEL, HOFF, NKEEP, REBT, T0, tag, dc.L, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_Peven_partialhankel_%s_off%d_reb%d_T0%d_%s_L%d_claude.png" % (GLABSEL.replace("+", "-"), HOFF, NKEEP, T0, tag, dc.L)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
