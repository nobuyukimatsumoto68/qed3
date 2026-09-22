#!/usr/bin/env python3
# sigma2_smear_gevp_v2_claude.py  [chunk 1: N_v-truncation (distillation-smearing) PS^2 point-operator matrix]
#   Expand the PS {sigma^2_00, O_2m, O_1m} point GEVP with a multi-smearing variational basis built by
#   TRUNCATING the distillation perambulator to N_v in SMEARS (no new solves -- block-slice the full peram).
#   Smearing to n modes = projector P_n = V_n V_n^dag (V_n = V[:,:n], modes ordered by Wilson eigenvalue).
#   Leg smeared n at sink t', m at source t:
#       A_nm(t',t) = V[t'][:,:n] ( tau(t',t)[:n,:m] - (1/2) I[:n,:m] delta_{t't} ) V[t][:,:m]^dag
#   (the GW contact 1/2 is subtracted in MODE space on the equal-time diagonal -> position contact -1/2 P_min(n,m);
#    at n=24=full, V V^dag = I, reduces to the existing -1/2 I).  PS legs = forward improved; PS==FS for this
#   4-point (per-loop S+Stilde = 2x improved), so NCYC = 2^{#cyc} and the n=24 sub-block MUST reproduce
#   fs_channels_v2.  Refs: multi-smearing Morningstar-Peardon hep-lat/9901004; distillation Peardon 0905.2160.
#
#   Run (validate): VALIDATE=1 python3 sigma2_smear_gevp_v2_claude.py
#   Run (build 12-op cache): OMP_NUM_THREADS=1 NPROC=16 DTMAX=24 python3 sigma2_smear_gevp_v2_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
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
SMEARS = [int(x) for x in os.environ.get("SMEARS", "4,8,16,24").split(",")]
NS = dc.NS
OPLAB = ["sigma^2_00", "O_2m", "O_1m"]
OFF = [(0, 0), (0, 0), (0, SPLIT)]                 # (sink,source) time offsets per base op
NCYC = np.array([2.0 ** len(c) for c in G.PERMS])  # per-loop S+Stilde = 2x (PS==FS)
NOP = len(OPLAB)
NSM = len(SMEARS)
NTOT = NOP * NSM                                   # 12 operators (op x smear)


def make_smear_blocks(k):
    # returns AblkS(ta,tb,n,m) (cached) + twin, nsite ; A_nm = V_n tau[:n,:m] V_m^dag - (1/2)P on equal time
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    U = [V[tsrc0 + a].T for a in range(twin)]      # (2Ns, Nv)
    nsite = U[0].shape[0] // NS
    cache = {}

    def AblkS(ta, tb, n, m):
        key = (ta, tb, n, m)
        if key not in cache:
            core = tau[ta, tb][:n, :m].astype(complex).copy()
            if ta == tb:
                kk = min(n, m)
                core[np.arange(kk), np.arange(kk)] -= 0.5
            A = U[ta][:, :n] @ core @ U[tb][:, :m].conj().T
            cache[key] = A.reshape(nsite, NS, nsite, NS)
        return cache[key]

    return AblkS, twin, nsite


def perm_contrib_smear(cycles, vtime, vsmear, bA, vspec, Pmap):
    # as G.perm_contrib_folded but each leg sliced by the smearings of the two vertices it connects
    arrs = []
    idxs = []
    for cyc in cycles:
        m = len(cyc)
        for a in range(m):
            va = cyc[a]
            vb = cyc[(a + 1) % m]
            arr = bA[(vtime[va], vtime[vb], vsmear[va], vsmear[vb])]
            if vspec[va][2]:
                arr = arr[:, Pmap, :, :, :]
            if vspec[vb][2]:
                arr = arr[:, :, :, Pmap, :]
            arrs.append(arr)
            idxs.append('z' + vspec[va][0] + G.SPIN[va] + vspec[vb][0] + G.SPIN[vb])
    for v in range(4):
        wv = vspec[v][1]
        if wv is not None:
            arrs.append(wv)
            idxs.append(vspec[v][0])
    sub = ','.join(idxs) + '->'
    path = G._PATH_CACHE.get(sub)
    if path is None:
        path = np.einsum_path(sub, *arrs, optimize='optimal')[0]
        G._PATH_CACHE[sub] = path
    val = np.einsum(sub, *arrs, optimize=path)
    return ((-1.0) ** len(cycles)) * val


def matrix_one_config(k, dual):
    AblkS, twin, nsite = make_smear_blocks(k)
    dualf = dual.astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    C = np.full((NTOT, NTOT, DTMAX), np.nan)
    omax = max(max(o) for o in OFF)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt + omax < twin and s + omax < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        # collect needed (c1,c2) time-offset pairs and (n,m) smear pairs
        toffs = set()
        for a in range(NOP):
            for b in range(NOP):
                vt = [dt + OFF[a][0], dt + OFF[a][1], OFF[b][0], OFF[b][1]]
                for va in range(4):
                    for vb in range(4):
                        toffs.add((vt[va], vt[vb]))
        bA = {}
        for (c1, c2) in toffs:
            for n in SMEARS:
                for m in SMEARS:
                    bA[(c1, c2, n, m)] = np.array([AblkS(s + c1, s + c2, n, m) for s in s0s])
        # NOTE: raw matrix is ASYMMETRIC (O_1m is time-split OFF=(0,SPLIT)); symmetrize at ENSEMBLE level
        # (as fs_channels_v2 does), NOT per-config -- so compute every (ia,ib).
        for a in range(NOP):
            for ni, nsm_a in enumerate(SMEARS):
                ia = a * NSM + ni
                for b in range(NOP):
                    for nj, nsm_b in enumerate(SMEARS):
                        ib = b * NSM + nj
                        vt = [dt + OFF[a][0], dt + OFF[a][1], OFF[b][0], OFF[b][1]]
                        vsm = [nsm_a, nsm_a, nsm_b, nsm_b]  # sink op a (0,1), source op b (2,3)
                        vspec = G.op_vspec(a, ('i', 'j'), dualf, wY) + G.op_vspec(b, ('k', 'l'), dualf, wY)
                        v = 0.0
                        for ip, cyc in enumerate(G.PERMS):
                            v += NCYC[ip] * perm_contrib_smear(cyc, vt, vsm, bA, vspec, Pmap)
                        C[ia, ib, dt] = v.real / len(s0s)
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
        # n=24 sub-block must reproduce fs_channels_v2 (PS==FS) ; compare on 1 config
        k = dc.KS[0]
        C = matrix_one_config(k, dual)
        i24 = SMEARS.index(24)
        idx24 = [a * NSM + i24 for a in range(NOP)]      # the 3 full-smear ops
        import fs_channels_v2_claude as fc
        KER, OFFfc = fc.kernels(dual)
        Cfc = fc.matrix_one_config(k, KER, OFFfc, dual)  # (3,3,DTMAX)
        print("# VALIDATE n=24 sub-block vs fs_channels_v2 (PS==FS), config k=%d" % k)
        print("#  dt |   max|C_smear[24] - fs_channels_v2|   (rel)")
        for dt in range(1, 8):
            sub = np.array([[C[idx24[a], idx24[b], dt] for b in range(NOP)] for a in range(NOP)])
            d = np.abs(sub - Cfc[:, :, dt])
            scale = np.abs(Cfc[:, :, dt]).max() + 1e-300
            print("#  %2d |   %.3e   (%.3e)" % (dt, d.max(), d.max() / scale))
        return

    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    CACHEDIR = "sigma2_smear_cache_claude"
    os.makedirs(CACHEDIR, exist_ok=True)
    cache = "%s/sigma2_smear_%s_%dcfg_sm%s_d%d_claude.npy" % (
        CACHEDIR, tag.replace(".", "p"), len(ks), "-".join(map(str, SMEARS)), SPLIT)
    if os.path.exists(cache):
        allC = np.load(cache)
        print("# loaded cache <- %s (%d cfg)" % (cache, allC.shape[0]))
        return
    print("# computing %d configs, %d workers, SMEARS=%s -> %d ops ..." % (len(ks), NPROC, SMEARS, NTOT))
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
