#!/usr/bin/env python3
# sigma2_point_test_claude.py  -- look at a SINGLE point-sigma^2 operator's two-point correlator.
#   O_point = sigma(x1) sigma(x2)  with x1,x2 two FIXED degree-5 vertices (NO spatial sum, NO A_x Y00) --
#   built via fs_gevp_point.op_vspec_point (one-hot weights).  PP flavor (both PS-furnished) by default.
#   Diagonal two-point <O_point(t) O_point(0)> (point sink AND point source), jackknife effmass.
#   Purpose: see what a raw unaveraged point op gives (Fin's warning: individually noisy) before the GEVP.
#   Run: ENS=<..> NVDIR=<..> LREF=<1|2> X1=0 X2=6 FLAV=PP BINSIZE=10 python3 sigma2_point_test_claude.py

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

X1 = int(os.environ.get("X1", "0"))
X2 = int(os.environ.get("X2", "6"))
FLAV = os.environ.get("FLAV", "PP")
DTMAX = int(os.environ.get("DTMAX", "20"))
NCFG = int(os.environ.get("NCFG", "0"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
FLAVMASK = {"PP": [0, 0], "FF": [1, 1], "FP": [1, 0]}


def flavfac(fsmask):
    out = np.zeros(len(G.PERMS))
    for ip, cycles in enumerate(G.PERMS):
        f = 1.0
        for cyc in cycles:
            f *= (1.0 + (-1.0) ** sum(fsmask[v] for v in cyc))
        out[ip] = f
    return out


def point_corr_one_config(k, x1, x2, ffac):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    Pmap = G.antipodal_map()
    vspec = G.op_vspec_point(x1, x2, ('i', 'j'), nsite) + G.op_vspec_point(x1, x2, ('k', 'l'), nsite)
    C = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        vt = [dt, dt, 0, 0]                                    # sink sigma^2 at dt, source at 0 (no time split)
        offs = set((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        base = np.array([G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap) for cyc in G.PERMS])
        C[dt] = (ffac @ base).real / len(s0s)
    return C


def main():
    tag = dc.ENS.split("nu0")[0]
    ffac = flavfac(FLAVMASK[FLAV] + FLAVMASK[FLAV])
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    print("# ENS=%s L=%d NVDIR=%s  point-sigma^2 %s at fixed sites (x1=%d,x2=%d)  %d cfg"
          % (tag, dc.L, os.environ["NVDIR"], FLAV, X1, X2, len(ks)))
    allC = np.array([point_corr_one_config(k, X1, X2, ffac) for k in ks])
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    cen = blk.mean(0)
    jk = np.array([np.delete(blk, i, 0).mean(0) for i in range(nb)])
    cen_err = np.sqrt((nb - 1) * np.mean((jk - jk.mean(0)) ** 2, 0))
    print("\n#  dt |    C(err)              S/N")
    for dt in range(DTMAX):
        sn = cen[dt] / cen_err[dt] if cen_err[dt] > 0 else 0.0
        print("#  %2d | % .4e(%.1e)  %6.2f" % (dt, cen[dt], cen_err[dt], sn))
    with np.errstate(all="ignore"):
        em = np.log(jk[:, :-1] / jk[:, 1:])
    em_c = em.mean(0)
    em_e = np.sqrt((nb - 1) * np.mean((em - em_c) ** 2, 0))
    print("\n#  dt | a_t m_eff(err)")
    for dt in range(DTMAX - 1):
        if np.isfinite(em_c[dt]) and np.isfinite(em_e[dt]):
            print("#  %2d | %7.4f(%.4f)" % (dt, em_c[dt], em_e[dt]))


if __name__ == "__main__":
    main()
