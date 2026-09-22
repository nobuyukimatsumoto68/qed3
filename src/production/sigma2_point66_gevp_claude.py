#!/usr/bin/env python3
# sigma2_point66_gevp_claude.py -- FULL point-sigma^2 basis: all C(12,2)=66 pairs of the degree-5 vertices.
#   O_{(x1,x2)} = sigma(x1) sigma(x2) at FIXED vertex sites (NO sum, NO A_x Y00).  Builds the 66x66 correlator
#   matrix per config, then a PLAIN (no-Hankel, per the stress-tensor agent) rebased GEVP.  Idea: the point
#   ops add the RELATIVE-POSITION variational directions the l=0-summed basis lacks -> resolve the two-meson.
#   By icosahedral symmetry the 66 pairs carry only ~3 distinct separations (nearest x30, next x30, antipodal
#   x6); the rebase is meant to handle that degeneracy (else fall back to symmetry-averaged shells).
#
#   Efficient contraction: perm_contrib_open keeps the 4 operator-vertex sites (i,j sink; k,l source) as FREE
#   output indices over the 12 vertices -> a 12^4 tensor T[i,j,k,l] per (config,dt); the 66x66 matrix is
#   C[(i<j),(k<l)] = T[i,j,k,l].  Validated against fs_gevp_point.perm_contrib_folded with one-hot weights.
#   Run: ENS=.. NVDIR=.. LREF=1|2 FLAV=PP NKEEP=4 REBT=4 T0=3 BINSIZE=10 VALIDATE=1 python3 sigma2_point66_gevp_claude.py

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

FLAV = os.environ.get("FLAV", "PP")
DTMAX = int(os.environ.get("DTMAX", "20"))
NCFG = int(os.environ.get("NCFG", "0"))
NPROC = int(os.environ.get("NPROC", "1"))
NKEEP = int(os.environ.get("NKEEP", "4"))
REBT = int(os.environ.get("REBT", "4"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
VALIDATE = int(os.environ.get("VALIDATE", "0"))
FLAVMASK = {"PP": [0, 0], "FF": [1, 1], "FP": [1, 0]}
NS = G.NS
OUT = "ijkl"


def flavfac(fsmask):
    out = np.zeros(len(G.PERMS))
    for ip, cycles in enumerate(G.PERMS):
        f = 1.0
        for cyc in cycles:
            f *= (1.0 + (-1.0) ** sum(fsmask[v] for v in cyc))
        out[ip] = f
    return out


_POPEN = {}


def perm_contrib_open(cycles, vt, bA_v):
    # like perm_contrib_folded but NO weights and the 4 operator sites (i,j,k,l) kept as FREE output over the
    # vertex sub-lattice; returns (z, VN,VN,VN,VN).  Site letter G.SITE[v], spin G.SPIN[v] per vertex.
    arrs = []
    idxs = []
    for cyc in cycles:
        m = len(cyc)
        for a in range(m):
            va = cyc[a]
            vb = cyc[(a + 1) % m]
            arrs.append(bA_v[(vt[va], vt[vb])])
            idxs.append('z' + G.SITE[va] + G.SPIN[va] + G.SITE[vb] + G.SPIN[vb])
    sub = ','.join(idxs) + '->z' + OUT
    path = _POPEN.get(sub)
    if path is None:
        path = np.einsum_path(sub, *arrs, optimize='optimal')[0]
        _POPEN[sub] = path
    val = np.einsum(sub, *arrs, optimize=path)
    return ((-1.0) ** len(cycles)) * val


def tensor_one_config(k, V, ffac):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    VN = len(V)
    T = np.full((DTMAX, VN, VN, VN, VN), np.nan)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        vt = [dt, dt, 0, 0]
        offs = set((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        bA_v = {o: np.array([AblkS(s + o[0], s + o[1])[np.ix_(V, range(NS), V, range(NS))] for s in s0s]) for o in offs}
        acc = None
        for ip, cyc in enumerate(G.PERMS):
            if ffac[ip] == 0.0:
                continue
            contr = ffac[ip] * perm_contrib_open(cyc, vt, bA_v)
            acc = contr if acc is None else acc + contr
        T[dt] = acc.real.mean(0)                         # average over s0 (z), keep (VN^4)
    return T


def pairs_matrix(T, V):
    VN = len(V)
    prs = [(i, j) for i in range(VN) for j in range(i + 1, VN)]      # 66 unordered pairs
    npr = len(prs)
    C = np.full((npr, npr, DTMAX), np.nan)
    for a, (i, j) in enumerate(prs):
        for b, (kk, ll) in enumerate(prs):
            C[a, b] = T[:, i, j, kk, ll]
    return C, prs


_WK = {}


def _init(V, ffac):
    _WK["V"] = V
    _WK["ffac"] = ffac


def _work(k):
    return tensor_one_config(k, _WK["V"], _WK["ffac"])


def main():
    tag = dc.ENS.split("nu0")[0]
    ffac = flavfac(FLAVMASK[FLAV] + FLAVMASK[FLAV])
    V = G.five_vertices()
    print("# ENS=%s L=%d NVDIR=%s  FLAV=%s  degree-5 vertices=%d -> C(%d,2)=%d point ops"
          % (tag, dc.L, os.environ["NVDIR"], FLAV, len(V), len(V), len(V) * (len(V) - 1) // 2))

    if VALIDATE:
        k = dc.KS[0]
        T = tensor_one_config(k, V, ffac)
        # compare one entry T[i,j,k,l] to a direct one-hot perm_contrib_folded
        i, j, kk, ll = V.index(V[0]), 1, 2, 3
        AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
        Pmap = G.antipodal_map()
        dt = 4
        s0s = np.array([s for s in range(twin) if s + dt < twin])
        vt = [dt, dt, 0, 0]
        offs = set((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        vsp = G.op_vspec_point(V[i], V[j], ('i', 'j'), nsite) + G.op_vspec_point(V[kk], V[ll], ('k', 'l'), nsite)
        base = np.array([G.perm_contrib_folded(cyc, vt, bAS, vsp, Pmap) for cyc in G.PERMS])
        ref = (ffac @ base).real / len(s0s)
        got = T[dt, i, j, kk, ll]
        print("# [VALIDATE] one-hot perm_contrib_folded = %.6e ; open-tensor = %.6e ; rel diff = %.2e"
              % (ref, got, abs(ref - got) / (abs(ref) + 1e-30)))

    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    if NPROC > 1:
        import multiprocessing as mp
        with mp.Pool(NPROC, initializer=_init, initargs=(V, ffac)) as pool:
            allT = np.array(pool.map(_work, ks))
    else:
        allT = np.array([tensor_one_config(k, V, ffac) for k in ks])
    ncfg = allT.shape[0]

    # 66x66 per config
    Cs = np.array([pairs_matrix(allT[c], V)[0] for c in range(ncfg)])   # (ncfg, 66, 66, DTMAX)
    _, prs = pairs_matrix(allT[0], V)
    Cs = 0.5 * (Cs + np.swapaxes(Cs, 1, 2))                             # symmetrize
    nb = ncfg // BINSIZE
    blk = np.array([Cs[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    # PLAIN (no Hankel) rebased GEVP, off=[0]
    import hankel_rebase_scan_claude as hs

    def gevp(Cmat, Vfix):
        Cts = np.transpose(Cmat, (2, 0, 1))
        Big = hs.hankel_off(Cts, [0])
        Vp = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
        return hs.rebased_effmass_fixed(Big, Vp, T0), Vp

    em_c, Vfix = gevp(blk.mean(0), None)
    ems = np.array([gevp(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]
    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)) + "  [66-op point basis, plain off=0]")
    for t in range(T0, min(tmax, 18)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(NKEEP))))


if __name__ == "__main__":
    main()
