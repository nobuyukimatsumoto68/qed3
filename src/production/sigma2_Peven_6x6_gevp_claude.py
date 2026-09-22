#!/usr/bin/env python3
# sigma2_Peven_6x6_gevp_claude.py  [dedicated parity-even P+ 6x6 GEVP + labeled figure]
#   Parity-even sector of the flavor x geometry basis: {PP, FF} x {sigma^2_00, O_2m, O_1m} = 6 ops
#   (flat indices flavor*3+geom, flavor 0=PP 1=FF -> OPS = 0,1,2,3,4,5).  Reads the FULL 9-op cache
#   built by sigma2_flavorgeom_full_v2_claude.py, selects the P+ block, and runs the PRODUCTION
#   block-Hankel + staged-rebase GEVP (Dt=OFFSETS, rebase NKEEP@REBT, metric point T0).  Writes the
#   labeled figure figs/sigma2_Peven_6x6_reb<NKEEP>_<tag>_claude.png (distinct from the full-9op fig).
#
#   T0 (metric point) is the timeslice whose C(T0) defines the GEVP metric; must satisfy
#   0 <= T0 < REBT and C(T0) well-conditioned (small T0; canonical 3, earliest 2).
#
#   Refs: block-Hankel/GPOF Aubin-Orginos 1010.0202 ; GEVP Luscher-Wolff ; distillation Peardon
#         0905.2160 ; flavor-cross furnishing rule ps_fs_flavor_cross_impl_plan_claude.md ;
#         sigma^2-F^2 mixing Chester-Pufu 1603.05582.
#   Run: OPS=0,1,2,3,4,5 OFFSETS=0,2,4 REBT=5 NKEEP=5 T0=2 BINSIZE=10 python3 sigma2_Peven_6x6_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs

SPLIT = int(os.environ.get("SPLIT", "1"))
T0 = int(os.environ.get("T0", "2"))
REBT = int(os.environ.get("REBT", "5"))
NKEEP = int(os.environ.get("NKEEP", "5"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,2,4").split(",")]
OPS = [int(x) for x in os.environ.get("OPS", "0,1,2,3,4,5").split(",")]

FLABS = ["PP", "FF", "FP"]
GEOLAB = ["s2", "O2m", "O1m"]
OPLAB = ["%s-%s" % (f, g) for f in FLABS for g in GEOLAB]    # flavor*3+geom
STATE_LAB = ["ground", "2-meson", "2-meson", "noisy", "exc 2-meson", "state 5"]


def hankel_reb(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    CACHEDIR = "sigma2_flavor_cache_claude"
    cands = ["%s/sigma2_flavorgeom_FULL_%s_%dcfg_d%d_claude.npy" % (CACHEDIR, tag.replace(".", "p"), n, SPLIT)
             for n in (400, 200, 100)]
    cache = next((c for c in cands if os.path.exists(c)), None)
    if cache is None:
        import glob
        hits = sorted(glob.glob("%s/sigma2_flavorgeom_FULL_%s_*cfg_d%d_claude.npy" % (CACHEDIR, tag.replace(".", "p"), SPLIT)))
        if not hits:
            print("# no FULL cache found; run sigma2_flavorgeom_full_v2_claude.py first")
            return
        cache = hits[-1]
    allC = np.load(cache)
    print("# loaded <- %s  shape %s" % (cache, allC.shape))

    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))    # ensemble symmetrize (O_1m time-split)
    allC = allC[:, OPS][:, :, OPS]
    ncfg = allC.shape[0]
    print("# ops: %s  (Hankel Dt=%s reb%d@%d T0=%d)" % ([OPLAB[s] for s in OPS], OFFSETS, NKEEP, REBT, T0))

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = hankel_reb(blk.mean(0), None)
    if nb < 2:
        em_err = np.zeros_like(em_c)                 # single config (free limit): exact, no jackknife
    else:
        ems = np.array([hankel_reb(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
        em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]

    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)) + "  [ref 0.46/0.62/0.92]")
    for t in range(T0, min(tmax, 22)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_err[t, n]) for n in range(NKEEP))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.axhline(0.644, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(tmax * 0.72, 0.652, r"$2m_{PS}=0.644$", fontsize=9, color="gray")
    ax.axhline(0.46, color="gray", ls=":", lw=1, alpha=0.45)
    ax.text(tmax * 0.72, 0.468, r"$0.46$ ground", fontsize=9, color="gray")
    ax.axhline(0.92, color="gray", ls=":", lw=1, alpha=0.45)
    ax.text(tmax * 0.62, 0.928, r"$0.92$ excited two-meson", fontsize=9, color="gray")
    cols = ["tab:green", "tab:red", "tab:orange", "gray", "tab:blue", "tab:purple"]
    mkr = ["o", "s", "s", "x", "D", "v"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.30)
        lab = "state%d (%s)" % (n, STATE_LAB[n] if n < len(STATE_LAB) else "state %d" % n)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 6], marker=mkr[n % 6], ms=5, lw=1.1,
                    capsize=2.5, label=lab)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.3, 1.1)
    ax.set_xlim(T0, min(tmax, int(os.environ.get("TMAXPLOT", "14"))))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"P+ 6x6 {PP,FF}$\times$geom reb%d@%d T0=%d  %s L1 %dcfg" % (NKEEP, REBT, T0, tag, ncfg), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_Peven_6x6_reb%d_T0%d_%s_claude.png" % (NKEEP, T0, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
