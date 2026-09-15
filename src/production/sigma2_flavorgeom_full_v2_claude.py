#!/usr/bin/env python3
# sigma2_flavorgeom_full_v2_claude.py  [chunk 2: FULL 9-op flavor x geometry GEVP]
#   Basis = {PP,FF,FP} (flavor) x {sigma^2_00, O_2m, O_1m} (geometry) = 9 operators.  Op index = flavor*3+geom.
#   Correlator C[(fa,ga),(fb,gb)] = FLAVFAC[fa,fb] @ base[ga,gb], where base[ga,gb][ip] = perm_contrib_folded
#   (forward improved, geometry (ga,gb)) and FLAVFAC[fa,fb][ip] = prod over cycles (1+(-1)^{sum fsmask}), fsmask
#   from the flavor pattern FLAV[fa] (vertices 0,1) + FLAV[fb] (vertices 2,3).  Production Hankel+rebase.
#   Combines the two variational axes (flavor furnishing + geometry profile).
#   Run: OFFSETS=0,2,4 REBT=4 NKEEP=3 T0=3 BINSIZE=10 NPROC=12 python3 sigma2_flavorgeom_full_v2_claude.py
#        OPS=... to select a sub-basis (flat indices flavor*3+geom).

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import glob
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import hankel_rebase_scan_claude as hs

SPLIT = int(os.environ.get("SPLIT", "1"))
DTMAX = int(os.environ.get("DTMAX", "24"))
NPROC = int(os.environ.get("NPROC", "1"))
NCFG = int(os.environ.get("NCFG", "0"))
T0 = int(os.environ.get("T0", "3"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,2,4").split(",")]
OPS = os.environ.get("OPS", ",".join(str(i) for i in range(9)))
OFF_GEOM = [(0, 0), (0, 0), (0, SPLIT)]
FLAV = {"PP": [0, 0], "FF": [1, 1], "FP": [1, 0]}
FLABS = ["PP", "FF", "FP"]
GEOLAB = ["s2", "O2m", "O1m"]
OPLAB = ["%s-%s" % (f, g) for f in FLABS for g in GEOLAB]    # flavor*3+geom


def flavfac(fsmask):
    out = np.zeros(len(G.PERMS))
    for ip, cycles in enumerate(G.PERMS):
        f = 1.0
        for cyc in cycles:
            f *= (1.0 + (-1.0) ** sum(fsmask[v] for v in cyc))
        out[ip] = f
    return out


FFAC = {(fa, fb): flavfac(FLAV[fa] + FLAV[fb]) for fa in FLABS for fb in FLABS}


def matrix_one_config(k, dual):
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(k)
    dualf = dual.astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    C = np.full((9, 9, DTMAX), np.nan)
    omax = max(max(o) for o in OFF_GEOM)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt + omax < twin and s + omax < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        offs = set()
        for ga in range(3):
            for gb in range(3):
                vt = [dt + OFF_GEOM[ga][0], dt + OFF_GEOM[ga][1], OFF_GEOM[gb][0], OFF_GEOM[gb][1]]
                offs.update((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        base = {}
        for ga in range(3):
            for gb in range(3):
                vt = [dt + OFF_GEOM[ga][0], dt + OFF_GEOM[ga][1], OFF_GEOM[gb][0], OFF_GEOM[gb][1]]
                vspec = G.op_vspec(ga, ('i', 'j'), dualf, wY) + G.op_vspec(gb, ('k', 'l'), dualf, wY)
                base[(ga, gb)] = np.array([G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap) for cyc in G.PERMS])
        for ifa, fa in enumerate(FLABS):
            for ga in range(3):
                ia = ifa * 3 + ga
                for ifb, fb in enumerate(FLABS):
                    for gb in range(3):
                        ib = ifb * 3 + gb
                        C[ia, ib, dt] = (FFAC[(fa, fb)] @ base[(ga, gb)]).real / len(s0s)
    return C


_WK = {}


def _init(dual):
    _WK["dual"] = dual


def _work(k):
    return matrix_one_config(k, _WK["dual"])


def hankel_reb(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    CACHEDIR = "sigma2_flavor_cache_claude"
    os.makedirs(CACHEDIR, exist_ok=True)
    cache = "%s/sigma2_flavorgeom_FULL_%s_%dcfg_d%d_claude.npy" % (CACHEDIR, tag.replace(".", "p"), len(ks), SPLIT)
    if os.path.exists(cache):
        allC = np.load(cache)
        print("# loaded <- %s" % cache)
    else:
        print("# building FULL 9-op flavor x geometry, %d cfg, %d workers ..." % (len(ks), NPROC))
        if NPROC > 1:
            import multiprocessing as mp
            with mp.Pool(NPROC, initializer=_init, initargs=(dual,)) as pool:
                allC = np.array(pool.map(_work, ks))
        else:
            allC = np.array([matrix_one_config(k, dual) for k in ks])
        np.save(cache, allC)
        print("# cached -> %s  shape %s" % (cache, allC.shape))

    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    sel = [int(x) for x in OPS.split(",")]
    allC = allC[:, sel][:, :, sel]
    ncfg = allC.shape[0]
    print("# ops: %s  (Hankel Dt=%s reb%d@%d T0=%d)" % ([OPLAB[s] for s in sel], OFFSETS, NKEEP, REBT, T0))
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = hankel_reb(blk.mean(0), None)
    ems = np.array([hankel_reb(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]
    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)) + "  [ref 0.46 / 0.62]")
    for t in range(T0, min(tmax, 22)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_err[t, n]) for n in range(NKEEP))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ax.axhline(0.644, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(tmax * 0.55, 0.652, r"$2m_{PS}=0.644$", fontsize=9, color="gray")
    ax.axhline(0.46, color="tab:green", ls=":", lw=1, alpha=0.5)
    cols = ["tab:green", "tab:red", "tab:blue", "tab:purple", "tab:orange"]
    mkr = ["o", "s", "^", "D", "v"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.25)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 5], marker=mkr[n % 5], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(tmax, 20))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FULL flavor$\times$geom 9-op (nop=%d) Hankel reb%d@%d  %s L1 %dcfg" % (len(sel), NKEEP, REBT, tag, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_flavorgeom_full_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
