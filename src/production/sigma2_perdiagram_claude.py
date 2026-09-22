#!/usr/bin/env python3
# sigma2_perdiagram_claude.py -- decompose <sigma^2_00(t) sigma^2_00(0)> (PP) into its 10 Wick diagrams and
#   group by fermion-loop topology, to test NM's DIAGRAM A hypothesis: the SAME-TIMESLICE cross-contraction
#   (both sink bilinears tied into ONE fermion loop = a nonlocal single-meson kernel psibar(x)K(x,y)psi(y))
#   carries the SINGLE-MESON tower (ground ~0.46 + the both-legs-n=1 excited E=4), while the FACTORIZED
#   diagrams (two independent meson loops) carry the genuine two-meson.
#   Groups: SINGLE-LOOP (len(cycles)==1, one 4-cycle threading all 4 vertices = diagram-A class) vs
#           MULTI-LOOP (>=2 cycles = two independent mesons + partials).  Effmass each + the full sum.
#   Run: ENS=.. NVDIR=.. LREF=1|2 FLAV=PP BINSIZE=10 python3 sigma2_perdiagram_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
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
BINSIZE = int(os.environ.get("BINSIZE", "10"))
FLAVMASK = {"PP": [0, 0], "FF": [1, 1], "FP": [1, 0]}
NPERM = len(G.PERMS)
NCYC = np.array([len(c) for c in G.PERMS])          # number of fermion loops per diagram


def flavfac(fsmask):
    out = np.zeros(NPERM)
    for ip, cycles in enumerate(G.PERMS):
        f = 1.0
        for cyc in cycles:
            f *= (1.0 + (-1.0) ** sum(fsmask[v] for v in cyc))
        out[ip] = f
    return out


def perdiag_one_config(k):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    vspec = G.op_vspec(0, ('i', 'j'), dualf, wY) + G.op_vspec(0, ('k', 'l'), dualf, wY)   # sigma^2_00
    B = np.full((NPERM, DTMAX), np.nan)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        vt = [dt, dt, 0, 0]
        offs = set((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        B[:, dt] = np.array([G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap).real for cyc in G.PERMS]) / len(s0s)
    return B


def effmass(blk):
    nb = blk.shape[0]
    cen = blk.mean(0)
    jk = np.array([np.delete(blk, i, 0).mean(0) for i in range(nb)])
    with np.errstate(all="ignore"):
        em = np.log(jk[:, :-1] / jk[:, 1:])
    em_c = em.mean(0)
    em_e = np.sqrt((nb - 1) * np.mean((em - em_c) ** 2, 0))
    return em_c, em_e


def main():
    tag = dc.ENS.split("nu0")[0]
    ffac = flavfac(FLAVMASK[FLAV] + FLAVMASK[FLAV])
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    print("# ENS=%s L=%d  FLAV=%s  %d perms; ncyc=%s  %d cfg" % (tag, dc.L, FLAV, NPERM, list(NCYC), len(ks)))
    allB = np.array([perdiag_one_config(k) for k in ks])       # (ncfg, NPERM, DTMAX)
    ncfg = allB.shape[0]
    nb = ncfg // BINSIZE

    # weighted per-diagram; group by loop count
    W = allB * ffac[None, :, None]                              # FLAVFAC-weighted
    single = NCYC == 1                                          # one 4-cycle = diagram-A class (single loop)
    grps = {"SINGLE-LOOP (diagram-A class)": single, "MULTI-LOOP (two-meson + partials)": ~single, "FULL sum": np.ones(NPERM, bool)}
    for name, mask in grps.items():
        C = W[:, mask, :].sum(1)                                # (ncfg, DTMAX)
        blk = np.array([C[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
        em_c, em_e = effmass(blk)
        cen = blk.mean(0)
        print("\n## %s   ( nperm=%d )" % (name, int(mask.sum())))
        print("#  dt |   C            a_t m_eff(err)")
        for dt in range(min(DTMAX - 1, 16)):
            m = "%7.4f(%.4f)" % (em_c[dt], em_e[dt]) if np.isfinite(em_c[dt]) and np.isfinite(em_e[dt]) else "  --"
            print("#  %2d | % .3e   %s" % (dt, cen[dt], m))


if __name__ == "__main__":
    main()
