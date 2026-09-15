#!/usr/bin/env python3
# sigma2_flavor_cross_v2_claude.py  [chunk 1: flavor-cross matrix {PP, FF, FP} at a point geometry]
#   Flavor two-meson basis {PP=sigma_PS^2, FF=sigma_FS^2, FP=sigma_FS sigma_PS}.  Diagonal PP=FF (PS=FS proven),
#   but the cross <FF PP> and FP genuinely differ (a loop threading FS & PS vertices breaks the GW no-op).
#   RULE (derived, see ps_fs_flavor_cross_impl_plan_claude.md): each closed loop contributes a flavor factor
#       (1 + (-1)^{n_FS}) , n_FS = # FS vertices in the loop,
#   times the forward-improved loop.  Implementation = fs_channels_v2 (forward improved AblkS, 10 diagrams via
#   PERMS cycles) with the uniform NCYC=2^{#cyc} REPLACED by
#       FLAVFAC[ip] = prod over cycles c of (1 + (-1)^{ sum_{v in c} fsmask[v] }).
#   fsmask over 4 vertices (0,1 sink ; 2,3 source) set by the flavor pair: PP=[0,0], FF=[1,1], FP=[1,0] (1=FS).
#   PP -> 2^{#cyc} (= fs_channels_v2) ; FF -> = PP (PS=FS) ; cross FF-PP -> E(C_S^2) vanishes.  GEOM op default
#   sigma^2_00 (op 0).  Run: VALIDATE=1 python3 sigma2_flavor_cross_v2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G

SPLIT = int(os.environ.get("SPLIT", "1"))
DTMAX = int(os.environ.get("DTMAX", "24"))
NPROC = int(os.environ.get("NPROC", "1"))
NCFG = int(os.environ.get("NCFG", "0"))
VALIDATE = int(os.environ.get("VALIDATE", "0"))
GEOM = int(os.environ.get("GEOM", "0"))            # geometry op: 0=sigma^2_00, 1=O_2m, 2=O_1m
OFF_GEOM = [(0, 0), (0, 0), (0, SPLIT)]            # (sink,source) time offsets per geometry op
FLAV = {"PP": [0, 0], "FF": [1, 1], "FP": [1, 0]}  # fsmask per operator's 2 vertices (1=FS)
FLABS = ["PP", "FF", "FP"]


def flavfac(fsmask):
    # per-perm flavor factor array: prod over cycles of (1 + (-1)^{sum fsmask over cycle})
    out = np.zeros(len(G.PERMS))
    for ip, cycles in enumerate(G.PERMS):
        f = 1.0
        for cyc in cycles:
            nfs = sum(fsmask[v] for v in cyc)
            f *= (1.0 + (-1.0) ** nfs)
        out[ip] = f
    return out


# precompute FLAVFAC for each flavor pair (sink a, source b)
FF_CACHE = {}
for a in FLABS:
    for b in FLABS:
        fsmask = FLAV[a] + FLAV[b]                  # vertices [0,1]=sink a, [2,3]=source b
        FF_CACHE[(a, b)] = flavfac(fsmask)


def matrix_one_config(k, dual):
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(k)   # forward improved AblkS (full Nv)
    dualf = dual.astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    nf = len(FLABS)
    C = np.full((nf, nf, DTMAX), np.nan)
    og = OFF_GEOM[GEOM]
    omax = max(og)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt + omax < twin and s + omax < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        vt = [dt + og[0], dt + og[1], og[0], og[1]]
        offs = set((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        vspec = G.op_vspec(GEOM, ('i', 'j'), dualf, wY) + G.op_vspec(GEOM, ('k', 'l'), dualf, wY)
        # forward-improved per-perm contractions (flavor-independent base), then weight by FLAVFAC
        base = np.array([G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap) for cyc in G.PERMS])  # (nperm,)
        for ia, a in enumerate(FLABS):
            for ib, b in enumerate(FLABS):
                C[ia, ib, dt] = (FF_CACHE[(a, b)] @ base).real / len(s0s)
    return C


_WK = {}


def _init(dual):
    _WK["dual"] = dual


def _work(k):
    return matrix_one_config(k, _WK["dual"])


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()

    if VALIDATE:
        k = dc.KS[0]
        C = matrix_one_config(k, dual)
        import fs_channels_v2_claude as fc
        KER, OFFfc = fc.kernels(dual)
        Cfc = fc.matrix_one_config(k, KER, OFFfc, dual)   # (3,3,DTMAX) point ops; [GEOM,GEOM] = this geometry
        print("# VALIDATE flavor cross (GEOM=%d, %s) config k=%d" % (GEOM, ["s2_00", "O_2m", "O_1m"][GEOM], k))
        print("#  dt |  PP-fc(reldiff)   FF-PP(reldiff, =0 PS=FS)   <FF PP>  <PP PP>  E-vanish?")
        for dt in range(1, 9):
            pp = C[0, 0, dt]
            ff = C[1, 1, dt]
            cross = C[1, 0, dt]
            fcv = Cfc[GEOM, GEOM, dt]
            r_ppfc = abs(pp - fcv) / (abs(fcv) + 1e-300)
            r_ffpp = abs(ff - pp) / (abs(pp) + 1e-300)
            print("#  %2d |  %.3e    %.3e               % .3e  % .3e"
                  % (dt, r_ppfc, r_ffpp, cross, pp))
        return

    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    CACHEDIR = "sigma2_flavor_cache_claude"
    os.makedirs(CACHEDIR, exist_ok=True)
    cache = "%s/sigma2_flavor_%s_%dcfg_g%d_d%d_claude.npy" % (CACHEDIR, tag.replace(".", "p"), len(ks), GEOM, SPLIT)
    if os.path.exists(cache):
        print("# cache exists <- %s" % cache)
        return
    print("# computing %d configs, %d workers, GEOM=%d flavors %s ..." % (len(ks), NPROC, GEOM, FLABS))
    if NPROC > 1:
        import multiprocessing as mp
        with mp.Pool(NPROC, initializer=_init, initargs=(dual,)) as pool:
            allC = np.array(pool.map(_work, ks))
    else:
        allC = np.array([matrix_one_config(k, dual) for k in ks])
    np.save(cache, allC)
    print("# cached -> %s  shape %s" % (cache, allC.shape))


if __name__ == "__main__":
    main()
