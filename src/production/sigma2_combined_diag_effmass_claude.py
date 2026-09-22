#!/usr/bin/env python3
# sigma2_combined_diag_effmass_claude.py -- DIAGONAL effective masses of the expanded combined GEVP basis
#   (3 shell single-meson ops {ell1/2,ell3/2,ell5/2} + 3 sigma^2 two-meson ops {s2,O2m,O1m}).
#   Computes ONLY the 6 diagonal correlators C_ii(t) (skips the expensive cross blob), then plain
#   log-ratio effmass with binned jackknife errors.  Diagnostic companion to sigma2_combined_gevp_claude.py.
#   Run: MODE_CONTACT=1 ENS=<L2> LREF=2 NVDIR=distill_Nv24_v2 WINDOWS=0-4,4-12,12-24 OPS2=0,1,2 \
#        BINSIZE=10 NCFG=200 python3 sigma2_combined_diag_effmass_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MODE_CONTACT", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import sigma2_mPS_gevp_claude as MG
import sigma2_combined_gevp_claude as CB

WINDOWS = CB.WINDOWS
OPS2 = CB.OPS2
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
NS = G.NS

OPNAMES = ["shell l1/2", "shell l3/2 (2,2)", "shell l5/2"] + ["s2", "O2m", "O1m"]


def one_config_diag(k, dualf, wY, Pmap):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    nsh = len(WINDOWS)
    nop = len(OPS2)
    ntot = nsh + nop
    Ps = [[CB.shell_projector(tau[a, a] - 0.5 * np.eye(tau.shape[-1]), lo, hi) for (lo, hi) in WINDOWS] for a in range(twin)]
    C = np.full((ntot, DTMAX), np.nan)
    for dt in range(DTMAX):
        s0s = np.array([s for s in range(twin) if s + dt < twin])
        if len(s0s) == 0:
            continue
        # shell diagonals (mode space)
        for a in range(nsh):
            acc = 0.0
            for s in s0s:
                t = s + dt
                acc += -np.trace(Ps[t][a] @ tau[t, s] @ Ps[s][a] @ tau[s, t]).real
            C[a, dt] = acc / len(s0s)
        # sigma^2 diagonals (PP 4-vertex)
        offs = {(0, 0), (dt, dt), (0, dt), (dt, 0)}
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        vt = [dt, dt, 0, 0]
        for ib, op in enumerate(OPS2):
            vspec = G.op_vspec(op, ('i', 'j'), dualf, wY) + G.op_vspec(op, ('k', 'l'), dualf, wY)
            base = np.array([G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap) for cyc in G.PERMS])
            C[nsh + ib, dt] = (MG.FFPP @ base).real / len(s0s)
    return C


def main():
    tag = dc.ENS.split("nu0")[0]
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    print("# ENS=%s L=%d  DIAGONAL effmass  windows=%s ops2=%s  %d cfg  MC=%s"
          % (tag, dc.L, WINDOWS, OPS2, len(ks), os.environ.get("MODE_CONTACT")))
    allC = np.array([one_config_diag(k, dualf, wY, Pmap) for k in ks])
    ncfg = allC.shape[0]
    nb = max(ncfg // BINSIZE, 1)
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    Cmean = blk.mean(0)
    # plain log-ratio effmass with jackknife over bins
    def effmass(Cm):
        return np.log(Cm[:, :-1] / Cm[:, 1:])
    em_c = effmass(Cmean)
    if nb < 2:
        em_e = np.zeros_like(em_c)
    else:
        ems = np.array([effmass(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
        em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    ntot = Cmean.shape[0]
    print("\n#  t | " + " | ".join("%-16s" % OPNAMES[i] for i in range(ntot)))
    for t in range(1, min(DTMAX - 1, 15)):
        row = " | ".join("%7.4f(%.4f)" % (em_c[i, t], em_e[i, t]) for i in range(ntot))
        print("#  %2d | %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    ts = np.arange(em_c.shape[1])
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.7)
    ax.text(13.5, m2ps + 0.01, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="k", ls=":", lw=1.0, alpha=0.6)
    ax.text(13.5, mps + 0.01, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="k")
    cols = ["firebrick", "tab:purple", "tab:blue", "tab:green", "tab:orange", "tab:brown"]
    mkr = ["o", "D", "s", "^", "v", "P"]
    for i in range(ntot):
        g = np.isfinite(em_c[i]) & np.isfinite(em_e[i]) & (em_e[i] < 0.4)
        ax.errorbar(ts[g], em_c[i, g], yerr=em_e[i, g], color=cols[i % 6], marker=mkr[i % 6], ms=5, lw=1.1,
                    capsize=2.5, label=OPNAMES[i])
    ax.set_ylim(0.2, 1.1)
    ax.set_xlim(1, 15)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"DIAGONAL effmass, combined basis  %s L%d %dcfg  win=%s ops2=%s"
                 % (tag, dc.L, ncfg, WINDOWS, OPS2), fontsize=9)
    ax.legend(fontsize=9, loc="upper right", ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_combined_diag_effmass_%s_L%d_claude.png" % (tag, dc.L)
    fig.savefig(out, dpi=140)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
