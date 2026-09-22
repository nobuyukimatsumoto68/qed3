#!/usr/bin/env python3
# sigma2_shell_Uwindow_claude.py -- shell single-meson operators with the kernel built from the DISTILLATION
#   EIGENVECTORS U (index windows), NOT from M = tau(t,t) - 1/2.
#   Kernel K_w(t) = sum_{i in w} u_i u_i^H ; the peram is already in the U basis, so the mode-space vertex is the
#   0/1 diagonal block E_w and
#     C_ab(t,s) = -Re sum_{i in a} sum_{j in b} tau(t,s)_ij tau(s,t)_ji .
#   NVDIR selects the basis: distill_Nv24_sym (evecs of Dtilde^H Dtilde + Dtilde Dtilde^H) or distill_Nv24.
#   WINDOWS = slicing pattern (index windows, ascending eigenvalue), independent of M.
#   Refs: distillation Peardon et al. arXiv:0905.2160.
#   Run: ENS=<L1> LREF=1 NVDIR=distill_Nv24_sym WINDOWS=0-4,4-12,12-24 NKEEP=3 REBT=3 T0=2 BINSIZE=10 NCFG=200 \
#        python3 sigma2_shell_Uwindow_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import h5py
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import hankel_rebase_scan_claude as hs

WINDOWS = [tuple(int(x) for x in w.split("-")) for w in os.environ.get("WINDOWS", "0-4,4-12,12-24").split(",")]
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
REBT = int(os.environ.get("REBT", "3"))
T0 = int(os.environ.get("T0", "2"))
NKEEP = int(os.environ.get("NKEEP", str(len(WINDOWS))))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))


def check_evals_ascending(k):
    fn = dc.PERAM_DIR + "peram.%d.h5" % k
    h = h5py.File(fn, "r")
    ev = np.array(h["evals"])
    h.close()
    assert np.all(np.diff(ev, axis=1) >= -1e-12), "evals not ascending: index windows would be wrong"


def one_config(k):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    nsh = len(WINDOWS)
    idx = [np.arange(lo, hi) for (lo, hi) in WINDOWS]
    C = np.full((nsh, nsh, DTMAX), np.nan)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if len(s0s) == 0:
            continue
        for a in range(nsh):
            for b in range(nsh):
                acc = 0.0
                for s in s0s:
                    t = s + dt
                    X = tau[t, s][np.ix_(idx[a], idx[b])]
                    Y = tau[s, t][np.ix_(idx[b], idx[a])]
                    acc += -np.sum(X * Y.T).real
                C[a, b, dt] = acc / len(s0s)
    return C


def gevp(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Vp = hs.staged_project(Cts, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Cts, Vp, T0), Vp


def diag_effmass(Cmat):
    d = np.array([Cmat[i, i] for i in range(Cmat.shape[0])])
    with np.errstate(all="ignore"):
        return np.log(d[:, :-1] / d[:, 1:])


def jk(fn, blk):
    nb = blk.shape[0]
    c = fn(blk.mean(0))
    if nb < 2:
        return c, np.zeros_like(c)
    sm = np.array([fn(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
    return c, np.sqrt((nb - 1) * np.mean((sm - sm.mean(0)) ** 2, axis=0))


def main():
    tag = dc.ENS.split("nu0")[0]
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    check_evals_ascending(ks[0])
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# ENS=%s L=%d NVDIR=%s  U-INDEX-WINDOW shell kernel (no M)  windows=%s  %d cfg  reb%d@%d T0=%d"
          % (tag, dc.L, dc.NVDIR, WINDOWS, len(ks), REBT, NKEEP, T0))
    allC = np.array([one_config(k) for k in ks])
    ncfg = allC.shape[0]
    nb = max(ncfg // BINSIZE, 1)
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    de_c, de_e = jk(diag_effmass, blk)
    Vfix = gevp(blk.mean(0), None)[1]
    ge_c, ge_e = jk(lambda Cm: gevp(Cm, Vfix)[0], blk)
    nsh = len(WINDOWS)
    print("\n# DIAGONAL effmass")
    print("#  t | " + " | ".join("win %-11s" % ("[%d,%d)" % w) for w in WINDOWS))
    for t in range(1, 15):
        print("#  %2d | %s" % (t, " | ".join("%7.4f(%.4f)" % (de_c[i, t], de_e[i, t]) for i in range(nsh))))
    print("\n# shell GEVP effmass")
    for t in range(T0, min(ge_c.shape[0], 15)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (ge_c[t, n], ge_e[t, n]) for n in range(ge_c.shape[1]))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = ["firebrick", "tab:blue", "tab:green", "tab:purple", "tab:orange"]
    mkr = ["o", "s", "^", "D", "v"]
    os.makedirs("figs", exist_ok=True)
    for kind in ("diag", "gevp"):
        fig, ax = plt.subplots(figsize=(9.2, 6.0))
        ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.7)
        ax.text(TMAXPLOT * 0.7, m2ps + 0.01, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
        ax.axhline(mps, color="k", ls=":", lw=1.0, alpha=0.6)
        ax.text(TMAXPLOT * 0.7, mps + 0.01, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="k")
        if kind == "diag":
            for i in range(nsh):
                ts = np.arange(de_c.shape[1])
                g = np.isfinite(de_c[i]) & np.isfinite(de_e[i]) & (de_e[i] < 0.3)
                ax.errorbar(ts[g], de_c[i, g], yerr=de_e[i, g], color=cols[i % 5], marker=mkr[i % 5], ms=5, lw=1.1,
                            capsize=2.5, label="window [%d,%d)" % WINDOWS[i])
        else:
            for n in range(ge_c.shape[1]):
                ts = np.arange(ge_c.shape[0])
                g = np.isfinite(ge_c[:, n]) & np.isfinite(ge_e[:, n]) & (ge_e[:, n] < 0.3)
                ax.errorbar(ts[g], ge_c[g, n], yerr=ge_e[g, n], color=cols[n % 5], marker=mkr[n % 5], ms=5, lw=1.1,
                            capsize=2.5, label="state %d" % n)
        ax.set_ylim(0.2, 1.1)
        ax.set_xlim(1, TMAXPLOT)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
        ax.set_title(r"U-index-window shell kernel (%s)  %s L%d %s %dcfg  win=%s"
                     % (kind, tag, dc.L, dc.NVDIR, ncfg, WINDOWS), fontsize=9)
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        out = "figs/sigma2_shell_Uwindow_%s_%s_L%d_%s_claude.png" % (kind, tag, dc.L, dc.NVDIR)
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print("# -> %s" % out)


if __name__ == "__main__":
    main()
