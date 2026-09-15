#!/usr/bin/env python3
# sigma2_flavorgeom_full_v2_nsrc2_claude.py  [nsrc=2 DATA variant of sigma2_flavorgeom_full_v2_claude.py]
#   Builds the FULL 9-op flavor x geometry correlator cache on the nsrc=2 perambulators
#   (NVDIR=distill_Nv24_v2, source windows tsrc_list=[0,64]).  The sigma^2 four-point is EQUAL-TIME within
#   a window, so each window w is contracted SEPARATELY (make_config_win builds AblkS from tau_w alone) and
#   ONLY the finished per-window correlator matrices C_w[9,9,dt] are averaged at fixed separation:
#       C(dt) = (1/n_win) sum_w C_w(dt).    tau is NEVER combined across windows.
#   For nsrc=1 files (v1) this reduces to the single-window build -> IDENTICAL to sigma2_flavorgeom_full_v2
#   (validation: run with NVDIR=distill_Nv24 and the cache matches v1 per config).
#   Cache -> sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_<tag>_<ncfg>cfg_nsrc2_d<SPLIT>_claude.npy.
#   Run: NVDIR=distill_Nv24_v2 NPROC=12 python3 sigma2_flavorgeom_full_v2_nsrc2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24_v2")
import sys
sys.path.insert(0, ".")
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


def matrix_one_window(AblkS, twin, dualf, wY, Pmap):
    # the full 9x9 flavor x geometry correlator for ONE source window (AblkS built from that window's tau)
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


def matrix_one_config(k, dual):
    # POST-CONTRACTION window average: contract the 9x9 matrix fully within each window, then mean over windows
    dualf = dual.astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    V, windows = dc.load_peram_windows(k)
    nw = len(windows)
    acc = None
    for w in range(nw):
        AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, w)
        Cw = matrix_one_window(AblkS, twin, dualf, wY, Pmap)
        acc = Cw if acc is None else acc + Cw
    return acc / nw


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
    cache = "%s/sigma2_flavorgeom_FULL_%s_%dcfg_nsrc2_d%d_claude.npy" % (CACHEDIR, tag.replace(".", "p"), len(ks), SPLIT)
    if os.path.exists(cache):
        allC = np.load(cache)
        print("# loaded <- %s" % cache)
    else:
        print("# building FULL 9-op flavor x geometry (nsrc2, NVDIR=%s), %d cfg, %d workers ..."
              % (os.environ["NVDIR"], len(ks), NPROC))
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


if __name__ == "__main__":
    main()
