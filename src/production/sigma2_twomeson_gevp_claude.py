#!/usr/bin/env python3
# sigma2_twomeson_gevp_claude.py -- two-meson-only GEVP {s2, O2m, O1m} (PP flavor, MODE_CONTACT) on a chosen
#   perambulator basis (NVDIR).  Lean: builds only the 3x3 sigma^2 x sigma^2 block per config (no shell blob),
#   multiprocessing over configs (NPROC x OMP=1), cache keyed by NVDIR.  Then per-mode Hankel + rebase + GEVP.
#   Two settings per run: (A) no Hankel NKEEP=3 ; (B) partial Hankel O2m,O1m [0,2], reb@4.
#   Refs: block-Hankel/GPOF Aubin-Orginos arXiv:1010.0202; distillation Peardon et al. arXiv:0905.2160.
#   Run: ENS=<L2> LREF=2 NVDIR=distill_Nv24_sym NCFG=200 python3 sigma2_twomeson_gevp_claude.py

import os
os.environ.setdefault("MODE_CONTACT", "1")
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
import sys
sys.path.insert(0, ".")
import numpy as np
from multiprocessing import Pool
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import sigma2_mPS_gevp_claude as MG
import hankel_rebase_scan_claude as hs

OPS2 = [0, 1, 2]
OPNAMES = ["s2", "O2m", "O1m"]
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
NPROC = int(os.environ.get("NPROC", "4"))
REBT = int(os.environ.get("REBT", "3"))
T0 = int(os.environ.get("T0", "2"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "16"))
SETTINGS = [("A_noHankel", [[0], [0], [0]], 3), ("B_partialHankel", [[0], [0, 2], [0, 2]], 4)]


def one_config(k):
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    nop = len(OPS2)
    C = np.full((nop, nop, DTMAX), np.nan)
    for dt in range(DTMAX):
        s0s = np.array([s for s in range(twin) if s + dt < twin])
        if len(s0s) == 0:
            continue
        offs = {(0, 0), (dt, dt), (0, dt), (dt, 0)}
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        vt = [dt, dt, 0, 0]
        for ia, opa in enumerate(OPS2):
            for ib, opb in enumerate(OPS2):
                vspec = G.op_vspec(opa, ('i', 'j'), dualf, wY) + G.op_vspec(opb, ('k', 'l'), dualf, wY)
                base = np.array([G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap) for cyc in G.PERMS])
                C[ia, ib, dt] = (MG.FFPP @ base).real / len(s0s)
    return C


def gevp(Cmat, offs_list, nkeep):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_permode(Cts, offs_list)
    Vp = hs.staged_project(Big, [(REBT, nkeep)], T0)
    return hs.rebased_effmass_fixed(Big, Vp, T0)


def main():
    tag = dc.ENS.split("nu0")[0]
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# ENS=%s L=%d NVDIR=%s  TWO-MESON-ONLY GEVP  %d cfg  MC=%s" % (tag, dc.L, dc.NVDIR, len(ks), os.environ["MODE_CONTACT"]))
    os.makedirs("cache_claude", exist_ok=True)
    cf = "cache_claude/twomeson_%s_L%d_%s_mc%s_n%d.npy" % (tag, dc.L, dc.NVDIR, os.environ["MODE_CONTACT"], len(ks))
    if os.path.exists(cf):
        print("# loading cache <- %s" % cf)
        allC = np.load(cf)
    else:
        with Pool(NPROC) as pool:
            allC = np.array(pool.map(one_config, ks, chunksize=4))
        np.save(cf, allC)
        print("# cached -> %s" % cf)
    ncfg = allC.shape[0]
    nb = max(ncfg // BINSIZE, 1)
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs("figs", exist_ok=True)
    for name, offs, nkeep in SETTINGS:
        em_c = gevp(blk.mean(0), offs, nkeep)
        if nb < 2:
            em_e = np.zeros_like(em_c)
        else:
            ems = np.array([gevp(np.delete(blk, i, 0).mean(0), offs, nkeep) for i in range(nb)])
            em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
        print("\n# %s  offsets=%s nkeep=%d" % (name, offs, nkeep))
        for t in range(T0, min(em_c.shape[0], 15)):
            print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(em_c.shape[1]))))
        fig, ax = plt.subplots(figsize=(9.2, 6.0))
        ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.7)
        ax.text(TMAXPLOT * 0.7, m2ps + 0.01, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
        ax.axhline(mps, color="firebrick", ls=":", lw=1.2, alpha=0.8)
        ax.text(TMAXPLOT * 0.7, mps + 0.01, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="firebrick")
        cols = ["firebrick", "tab:purple", "tab:blue", "tab:green"]
        mkr = ["o", "D", "s", "^"]
        ts = np.arange(em_c.shape[0])
        for n in range(em_c.shape[1]):
            g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
            ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 4], marker=mkr[n % 4], ms=5, lw=1.1,
                        capsize=2.5, label="state %d" % n)
        ax.set_ylim(0.2, 1.3)
        ax.set_xlim(T0, TMAXPLOT)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
        ax.set_title("two-meson-only GEVP %s  offs=%s reb%d@%d T0=%d  %s L%d %s %dcfg"
                     % (name, offs, REBT, nkeep, T0, tag, dc.L, dc.NVDIR, ncfg), fontsize=8)
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        out = "figs/sigma2_twomeson_gevp_%s_%s_L%d_%s_claude.png" % (name, tag, dc.L, dc.NVDIR)
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print("# -> %s" % out)


if __name__ == "__main__":
    main()
