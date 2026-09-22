#!/usr/bin/env python3
# sigma2_Peven_partialhankel_fit_claude.py -- ONE partial-Hankel combination (s2 offset-set S2SET, O1m
#   offset-set O1MSET) + a CORRELATED constant (plateau) fit of each GEVP effmass over t in [TFIT_LO,TFIT_HI].
#   Fit = GLS with the jackknife covariance of the effmass points (u^T Cinv y / u^T Cinv u); stat error from
#   (u^T Cinv u)^{-1/2}; chi2/dof reported.  Uncorrelated (diagonal) fit also printed for cross-check.
#   Basis: P+ 6 base ops (off 0) + s2 replicas at nonzero offsets of S2SET + O1m replicas of O1MSET (both PP,FF).
#   Run: ENS=.. NVDIR=.. LREF=2 S2SET=0-2 O1MSET=0-2-4 TFIT_LO=5 TFIT_HI=9 REBT=4 NKEEP=3 T0=3 BINSIZE=10 \
#        python3 sigma2_Peven_partialhankel_fit_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24_v2")
import sys
sys.path.insert(0, ".")
import glob
import re
import numpy as np
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs

REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "3"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
SPLIT = int(os.environ.get("SPLIT", "1"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "16"))
TFIT_LO = int(os.environ.get("TFIT_LO", "5"))
TFIT_HI = int(os.environ.get("TFIT_HI", "9"))

IDX_S2 = [0, 3]                                             # PP-s2, FF-s2
IDX_O1M = [2, 5]                                            # PP-O1m, FF-O1m
BASE = [0, 1, 2, 3, 4, 5]
S2SET = [int(x) for x in os.environ.get("S2SET", "0-2").split("-")]
O1MSET = [int(x) for x in os.environ.get("O1MSET", "0-2-4").split("-")]


def build_basis():
    baseidx = list(BASE)
    off = [0] * 6
    for o in S2SET:
        if o == 0:
            continue
        baseidx += IDX_S2
        off += [o, o]
    for o in O1MSET:
        if o == 0:
            continue
        baseidx += IDX_O1M
        off += [o, o]
    return baseidx, off


def build_big(Cmat, baseidx, off):
    n = len(baseidx)
    omax = max(off)
    DT = Cmat.shape[2]
    tmax = DT - 2 * omax
    Big = np.full((tmax, n, n), np.nan)
    for t in range(tmax):
        for a in range(n):
            for b in range(n):
                Big[t, a, b] = Cmat[baseidx[a], baseidx[b], t + off[a] + off[b]]
    return Big


def gevp(Cmat, baseidx, off, Vfix):
    Big = build_big(Cmat, baseidx, off)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def corr_const_fit(ems, em_c, lo, hi, n):
    # ems: (nb, tmax, NKEEP) jackknife effmass samples ; correlated constant fit for state n over t in [lo,hi].
    nb = ems.shape[0]
    ts = np.arange(lo, hi + 1)
    y = em_c[ts, n]
    S = ems[:, ts, n]                                       # (nb, npts)
    ybar = S.mean(0)
    C = (nb - 1) * np.mean((S - ybar)[:, :, None] * (S - ybar)[:, None, :], axis=0)   # jk cov of the mean
    u = np.ones(len(ts))
    Cinv = np.linalg.inv(C)
    denom = u @ Cinv @ u
    p_corr = (u @ Cinv @ y) / denom
    err_corr = np.sqrt(1.0 / denom)
    r = y - p_corr * u
    chi2 = r @ Cinv @ r
    dof = len(ts) - 1
    # uncorrelated (diagonal) cross-check
    d = np.diag(C)
    w = 1.0 / d
    p_diag = np.sum(w * y) / np.sum(w)
    err_diag = np.sqrt(1.0 / np.sum(w))
    chi2_diag = np.sum((y - p_diag) ** 2 * w)
    return p_corr, err_corr, chi2 / dof, p_diag, err_diag, chi2_diag / dof, dof


def main():
    tag = dc.ENS.split("nu0")[0]
    CACHEDIR = "sigma2_flavor_cache_claude"
    MCTAG = "_mc1" if int(os.environ.get("MODE_CONTACT", "0")) else ""   # mode-space contact fix -> _mc1 cache
    hits = glob.glob("%s/sigma2_flavorgeom_FULL_%s_L%d_*cfg_nsrc2_d%d%s_claude.npy" % (CACHEDIR, tag.replace(".", "p"), dc.L, SPLIT, MCTAG))
    if not hits:
        hits = glob.glob("%s/sigma2_flavorgeom_FULL_%s_*cfg_d%d%s_claude.npy" % (CACHEDIR, tag.replace(".", "p"), SPLIT, MCTAG))
    cache = max(hits, key=lambda p: int(re.search(r"_(\d+)cfg", p).group(1)))
    allC = np.load(cache)
    print("# loaded <- %s  shape %s" % (cache, allC.shape))
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    baseidx, off = build_basis()
    print("# combo: s2=%s  O1m=%s  ->  %d ops ; reb%d@%d T0=%d ; fit t=[%d,%d]"
          % ("-".join(map(str, S2SET)), "-".join(map(str, O1MSET)), len(baseidx), NKEEP, REBT, T0, TFIT_LO, TFIT_HI))
    em_c, Vfix = gevp(blk.mean(0), baseidx, off, None)
    ems = np.array([gevp(np.delete(blk, i, 0).mean(0), baseidx, off, Vfix)[0] for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]

    fits = []
    print("\n# state | plateau_corr(err)  chi2/dof | plateau_diag(err) chi2/dof")
    for n in range(NKEEP):
        pc, ec, x2c, pd, ed, x2d, dof = corr_const_fit(ems, em_c, TFIT_LO, TFIT_HI, n)
        fits.append((pc, ec, x2c))
        print("#   %d   | %7.4f(%.4f)  %6.2f  | %7.4f(%.4f)  %6.2f   (dof=%d)"
              % (n, pc, ec, x2c, pd, ed, x2d, dof))

    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)))
    for t in range(T0, min(tmax, TMAXPLOT + 1)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(NKEEP))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m2ps = float(os.environ.get("M2PS", "0.7054"))
    cols = ["tab:green", "tab:red", "tab:blue"]
    mkr = ["o", "s", "^"]
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(TMAXPLOT * 0.62, m2ps + 0.008, r"$2m_{PS}=%.3f$ (L%d)" % (m2ps, dc.L), fontsize=8, color="gray")
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n], marker=mkr[n], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
        pc, ec, x2c = fits[n]
        ax.fill_between([TFIT_LO, TFIT_HI], pc - ec, pc + ec, color=cols[n], alpha=0.22)
        ax.plot([TFIT_LO, TFIT_HI], [pc, pc], color=cols[n], lw=1.6)
        ax.text(TFIT_HI + 0.15, pc, r"$%.4f(%.4f)$ [$\chi^2/\nu{=}%.1f$]" % (pc, ec, x2c), fontsize=8, color=cols[n], va="center")
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.axvspan(TFIT_LO, TFIT_HI, color="gray", alpha=0.06)
    ax.set_ylim(0.25, 1.05)
    ax.set_xlim(T0, TMAXPLOT)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"P+ 6op + Hankel  s2:%s O1m:%s  reb%d@%d T0=%d  fit[%d,%d]  %s L%d %dcfg"
                 % ("-".join(map(str, S2SET)), "-".join(map(str, O1MSET)), NKEEP, REBT, T0, TFIT_LO, TFIT_HI, tag, dc.L, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_Peven_partialhankel_fit_s2%s_o1m%s_%s_L%d_reb%d_T0%d_fit%d%d_claude.png" % (
        "".join(map(str, S2SET)), "".join(map(str, O1MSET)), tag, dc.L, NKEEP, T0, TFIT_LO, TFIT_HI)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
