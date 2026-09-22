#!/usr/bin/env python3
# t00_sigma2_cross_interacting_claude.py
#   INTERACTING-L1 variant of t00_sigma2_cross_claude.py (the {1,1,1,1} thread's free triangle).
#   QUESTION (NM discriminator): the free-L1 <T_00 sigma^2_P+> = 0 (machine precision, complete basis) --
#   is it a GENUINE selection rule (config-independent Wick identity, survives interactions) or the FREE
#   off-site-vs-on-site orthogonality (O_H perp on-site sigma densities, which BREAKS under interactions)?
#     - genuine  => cross ~ MACHINE-ZERO PER CONFIG (Wick/trace identity, gauge-independent)
#     - free-only => cross ~ O(sigma^2 scale) once the sink W carries the config spatial link phases
#   Same operator + same triangle as the free script; here the sink W = ti.build_W_gauge(sp) uses each
#   config's spatial link phases (ckpoint_lat.k), source Phi = V^dag diag(w00) V per config, then jackknife.
#   L1 Nv=24 = COMPLETE basis for ANY gauge config, so tau^dag = delta - tau stays exact -> clean test.
#   Run: ENS=Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000 NVDIR=distill_Nv24 \
#        LREF=1 python3 t00_sigma2_cross_interacting_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import geom_hopping_claude as gh
import t00_stress_ham_interacting_claude as ti

DTMAX = int(os.environ.get("DTMAX", "16"))
KMIN = int(os.environ.get("KMIN", "20"))
NCFG = int(os.environ.get("NCFG", "0"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
LREF = int(os.environ.get("LREF", "1"))
GEOM = os.environ.get("GEOM", "../../geometry/data/")
CONFIGDIR = os.environ.get("CONFIGDIR", dc.ENS)


def triangle(PWa_list, Phi, tau, leg_sink, twin, dtmax):
    # identical Wick structure to t00_sigma2_cross_claude.triangle: 4 classes (mult 2,2,1,1),
    # sink O_H at t, source sigma^2 at s (s<t), both sink terms tie with plain tau (leg_sink).
    Iv = np.eye(tau.shape[-1])
    tt = [tau[a, a] - 0.5 * Iv for a in range(twin)]          # source equal-time contact tau(s,s)-1/2
    out = np.zeros((4, dtmax))
    for dt in range(dtmax):
        ns = twin - dt
        acc = np.zeros(4)
        for s in range(ns):
            t = s + dt
            PW = PWa_list[t]
            Ph = Phi[s]
            tss = tt[s]
            Ls_ts = leg_sink[t, s]
            Ls_st = leg_sink[s, t]
            oaloop = np.trace(PW @ Ls_ts @ Ph @ Ls_st).real
            dS = np.trace(Ph @ tss).real
            dpS = np.trace(Ph @ tss @ Ph @ tss).real
            dW = np.trace(PW @ leg_sink[t, t]).real
            tri = -np.trace(PW @ Ls_ts @ Ph @ tss @ Ph @ Ls_st).real
            semi = oaloop * dS
            sinktad = dW * dpS
            disc = -dW * dS * dS
            acc += np.array([tri, semi, sinktad, disc])
        out[:, dt] = acc / ns
    return out


def corr_one(k, om, alpha, nns, nsite, tab, w00):
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    sp = ti.read_config_sp("%s/ckpoint_lat.%d" % (CONFIGDIR, k), len(tab) // 2)
    dtmax = min(DTMAX, twin)
    PhiW = []
    PhiWH = []
    Phi = []
    for a in range(twin):
        t = tsrc0 + a
        Vt = V[t].T
        W = ti.build_W_gauge(om, alpha, nns, nsite, tab, sp[t])       # gauge-dressed sink hop (r=0)
        PhiW.append(Vt.conj().T @ W @ Vt)
        PhiWH.append(Vt.conj().T @ W.conj().T @ Vt)
        Phi.append(Vt.conj().T @ (w00[:, None] * Vt))                 # sigma_PS source vertex
    triW = triangle(PhiW, Phi, tau, tau, twin, dtmax)
    triWH = triangle(PhiWH, Phi, tau, tau, twin, dtmax)
    mult = np.array([2.0, 2.0, 1.0, 1.0])
    cross = (mult[:, None] * triW).sum(0) + (mult[:, None] * triWH).sum(0)
    conn = 2.0 * triW[0] + 2.0 * triWH[0]
    # sigma 2pt norm reference (single-meson scale) + <T_00 T_00> sanity
    s2 = np.zeros(dtmax)
    chh = np.zeros(dtmax)
    for dt in range(dtmax):
        ns = twin - dt
        a2 = 0.0
        ah = 0.0
        for s in range(ns):
            t = s + dt
            a2 += -np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t]).real
            ah += -np.trace(PhiW[t] @ tau[t, s] @ PhiW[s] @ tau[s, t]).real
            ah += -np.trace(PhiWH[t] @ tau[t, s] @ PhiWH[s] @ tau[s, t]).real
        s2[dt] = a2 / ns
        chh[dt] = ah / ns
    return cross, conn, s2, chh


def main():
    om, alpha, nns, nsite = gh.build(GEOM, LREF)
    tab = ti.load_link_table("primal_links_n%d_claude.dat" % LREF)
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    ks = [k for k in dc.KS if k >= KMIN]
    if NCFG:
        ks = ks[:NCFG]
    print("# ENS=%s NVDIR=%s nsite=%d  ncfg=%d KMIN=%d  (interacting T_00 x sigma^2 cross)"
          % (dc.ENS, dc.NVDIR, nsite, len(ks), KMIN))
    allc = []
    allconn = []
    alls2 = []
    allh = []
    for i, k in enumerate(ks):
        cross, conn, s2, chh = corr_one(k, om, alpha, nns, nsite, tab, w00)
        allc.append(cross)
        allconn.append(conn)
        alls2.append(s2)
        allh.append(chh)
        if i == 0:
            r = np.abs(cross) / np.abs(s2)
            print("\n# ---- FIRST CONFIG k=%d (per-config decisive test) ----" % k)
            print("#  dt |    cross(k)      C_S(sigma2pt,k)   |cross|/|C_S|")
            for dt in range(1, min(12, cross.shape[0])):
                print("#  %2d | %13.5e   %13.5e   %10.3e" % (dt, cross[dt], s2[dt], r[dt]))
            mx = np.nanmax(r[1:11])
            print("#  => max |cross|/|C_S| over dt[1,10] = %.3e  (%s)"
                  % (mx, "MACHINE-ZERO -> genuine selection rule" if mx < 1e-8
                     else "O(1) -> NONZERO, free-only orthogonality broken by interactions"))
        if (i + 1) % 50 == 0:
            print("#   ... %d/%d" % (i + 1, len(ks)))
    allc = np.array(allc)
    allconn = np.array(allconn)
    alls2 = np.array(alls2)
    allh = np.array(allh)
    ncfg = allc.shape[0]
    dtmax = allc.shape[1]

    nb = max(ncfg // BINSIZE, 1)
    def jk(arr):
        blk = np.array([arr[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
        mean = blk.mean(0)
        jm = np.array([np.delete(blk, i, 0).mean(0) for i in range(nb)])
        err = np.sqrt((nb - 1) * np.mean((jm - jm.mean(0)) ** 2, axis=0)) if nb > 1 else np.zeros_like(mean)
        return mean, err
    cm, ce = jk(allc)
    s2m, s2e = jk(alls2)
    hm, he = jk(allh)

    print("\n# ==== JACKKNIFE over %d cfg (bin%d, nb=%d) ====" % (ncfg, BINSIZE, nb))
    print("# <T_00 T_00> sanity (interacting L1 T_00 ~ 0.459): C_HH(dt) plateau in effmass region")
    print("#  dt |   cross(mean +/- err)      |cross|/|C_S|(mean)   |  C_HH")
    for dt in range(1, dtmax):
        ratio = abs(cm[dt]) / abs(s2m[dt]) if s2m[dt] != 0 else np.nan
        nsig = abs(cm[dt]) / ce[dt] if ce[dt] > 0 else np.nan
        print("#  %2d | %12.4e +/- %10.3e  (%.1f sig)  %10.3e   | %11.4e"
              % (dt, cm[dt], ce[dt], nsig, ratio, hm[dt]))

    print("\n# VERDICT: cross consistent with ZERO within errors AND |cross|/|C_S| ~ machine eps => genuine;")
    print("#          cross many-sigma nonzero with |cross|/|C_S| = O(0.01-1) => free-only (interaction-broken).")


if __name__ == "__main__":
    main()
