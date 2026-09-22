#!/usr/bin/env python3
# sigma2_mPS_gevp_claude.py -- add the single-meson operator sigma_00 (ell=0 zero-momentum scalar bilinear) to
#   the P+ sigma^2 GEVP basis, so the GEVP assigns the (truncation-)leaked m_PS to its OWN state via sigma_00 and
#   leaves the two-meson to the sigma^2 operators.  Study the sigma^2 <-> m_PS channel (NM).
#
#   Basis (default): {sigma_00 (1 bilinear, 2 fermions), sigma^2_00 (PP, s2; 2 bilinears, 4 fermions)}.  2x2.
#   Correlators, ALL connected, ALL from the same peram + same wY weight + same loop-sign convention:
#     C_11 = <sigma_00(t) sigma_00(0)>       : 1 fermion loop, 2 vertices.
#     C_12 = <sigma_00(t) sigma^2_00(0)>     : 1 fermion loop, 3 vertices (the triangle).  C_21 = symmetrize.
#     C_22 = <sigma^2_00(t) sigma^2_00(0)>   : bridging perms (4 vertices), flavfac PP.
#   No sigma3 assumptions.  tau ONLY (never tau_gw; FS==PS collapse).  Contact tt = tau - 1/2 on equal time.
#   VALIDATE=1 checks C_22 here == flavorgeom cache C[0,0] (PP-s2) up to the machinery.
#   Run: ENS=.. LREF=2 NVDIR=distill_Nv24_v2 T0=2 REBT=4 NKEEP=2 BINSIZE=10 VALIDATE=1 python3 sigma2_mPS_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24_v2")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import hankel_rebase_scan_claude as hs

T0 = int(os.environ.get("T0", "2"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "2"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
DTMAX = int(os.environ.get("DTMAX", "20"))
NCFG = int(os.environ.get("NCFG", "0"))
VALIDATE = int(os.environ.get("VALIDATE", "0"))
NS = G.NS


def perm_contrib_n(cycles, vtime, bA, vspec, Pmap):
    # general n-vertex version of G.perm_contrib_folded (that one hardcodes 4 vertices).
    n = len(vspec)
    spin = [chr(ord('A') + i) for i in range(n)]
    arrs = []
    idxs = []
    for cyc in cycles:
        m = len(cyc)
        for a in range(m):
            va = cyc[a]
            vb = cyc[(a + 1) % m]
            arr = bA[(vtime[va], vtime[vb])]
            if vspec[va][2]:
                arr = arr[:, Pmap, :, :, :]
            if vspec[vb][2]:
                arr = arr[:, :, :, Pmap, :]
            arrs.append(arr)
            idxs.append('z' + vspec[va][0] + spin[va] + vspec[vb][0] + spin[vb])
    for v in range(n):
        wv = vspec[v][1]
        if wv is not None:
            arrs.append(wv)
            idxs.append(vspec[v][0])
    sub = ','.join(idxs) + '->'
    val = np.einsum(sub, *arrs, optimize='optimal')
    return ((-1.0) ** len(cycles)) * val


def flavfac_PP():
    out = np.zeros(len(G.PERMS))
    for ip, cyc in enumerate(G.PERMS):
        f = 1.0
        for c in cyc:
            f *= 2.0                     # all-PS: (1 + (-1)^0) = 2 per cycle
        out[ip] = f
    return out


FFPP = flavfac_PP()


def one_config(k, dualf, wY, Pmap):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    C = np.full((2, 2, DTMAX), np.nan)
    # vspecs (same wY convention as op_vspec).  single meson: 1 vertex; sigma^2_00: 2 vertices (s2 geom).
    vsp_ss = [(('i', wY, False)), (('k', wY, False))]                       # <sigma_00 sigma_00> : sink i, source k
    vsp_12 = [('i', wY, False)] + G.op_vspec(0, ('k', 'l'), dualf, wY)      # sink sigma_00 (i) + source sigma^2 (k,l)
    vsp_22 = G.op_vspec(0, ('i', 'j'), dualf, wY) + G.op_vspec(0, ('k', 'l'), dualf, wY)
    CYC_2 = [[0, 1]]
    CYC_3 = [[0, 1, 2], [0, 2, 1]]
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        # blocks needed: equal-time (s,s), (t,s), (s,t) with t=s+dt
        offs = {(0, 0), (dt, dt), (0, dt), (dt, 0)}
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        # C_11 : 2-vertex loop, vtime=[dt,0]
        C[0, 0, dt] = perm_contrib_n(CYC_2, [dt, 0], bAS, vsp_ss, Pmap).real / len(s0s)
        # C_12 : 3-vertex triangle, sink sigma_00 @ dt, source sigma^2 @ 0
        C[0, 1, dt] = perm_contrib_n(CYC_3, [dt, 0, 0], bAS, vsp_12, Pmap).real / len(s0s)
        # C_22 : 4-vertex, bridging perms, flavfac PP, vtime=[dt,dt,0,0]
        vt = [dt, dt, 0, 0]
        base = np.array([G.perm_contrib_folded(cyc, vt, bAS, vsp_22, Pmap) for cyc in G.PERMS])
        C[1, 1, dt] = (FFPP @ base).real / len(s0s)
    C[1, 0] = C[0, 1]                         # symmetric (real, translation-averaged)
    return C


def gevp(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, [0])
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    print("# ENS=%s L=%d NVDIR=%s  {sigma_00, sigma^2_00} 2x2  %d cfg  reb%d@%d T0=%d"
          % (tag, dc.L, os.environ["NVDIR"], len(ks), NKEEP, REBT, T0))
    allC = np.array([one_config(k, dualf, wY, Pmap) for k in ks])
    ncfg = allC.shape[0]

    if VALIDATE:
        import glob
        import re
        hits = glob.glob("sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_%s_L%d_*cfg_nsrc2_d1_claude.npy" % (tag.replace(".", "p"), dc.L))
        if not hits:
            hits = [h for h in glob.glob("sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_%s_*cfg_nsrc2_d1_claude.npy" % tag.replace(".", "p")) if "_L2_" not in h and "_L4_" not in h]
        if hits:
            cache = max(hits, key=lambda p: int(re.search(r"_(\d+)cfg", p).group(1)))
            cc = np.load(cache)
            ncmp = min(cc.shape[0], ncfg)
            r = allC[:ncmp, 1, 1, 2:8].mean(0) / cc[:ncmp, 0, 0, 2:8].mean(0)
            print("# [VALIDATE] C_22(here)/cache C[0,0] (PP-s2) dt=2..7 = %s  (should be a CONSTANT if consistent)"
                  % np.array2string(r, precision=4))

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = gevp(blk.mean(0), None)
    ems = np.array([gevp(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# m_PS=%.4f  2m_PS=%.4f" % (mps, m2ps))
    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)))
    for t in range(T0, min(tmax, 16)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(NKEEP))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(tmax * 0.62, m2ps + 0.008, r"$2m_{PS}=%.4f$ (L%d)" % (m2ps, dc.L), fontsize=8, color="gray")
    ax.axhline(mps, color="dimgray", ls=":", lw=1.1, alpha=0.6)
    ax.text(tmax * 0.62, mps + 0.008, r"$m_{PS}=%.4f$ (L%d)" % (mps, dc.L), fontsize=8, color="dimgray")
    cols = ["tab:green", "tab:red", "tab:blue"]
    mkr = ["o", "s", "^"]
    lab = ["state 0", "state 1", "state 2"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 3], marker=mkr[n % 3], ms=5, lw=1.1,
                    capsize=2.5, label=lab[n])
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.25, 1.1)
    ax.set_xlim(T0, min(tmax, int(os.environ.get("TMAXPLOT", "14"))))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"{$\sigma_{00}$, $\sigma^2_{00}$} 2x2 GEVP  reb%d@%d T0=%d  %s L%d %dcfg" % (NKEEP, REBT, T0, tag, dc.L, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_mPS_gevp_reb%d_T0%d_%s_L%d_claude.png" % (NKEEP, T0, tag, dc.L)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
