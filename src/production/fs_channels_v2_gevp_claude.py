#!/usr/bin/env python3
# fs_channels_v2_gevp_claude.py  [3x3 GEVP on the CORRECTED FS channel matrix {sigma^2_00, O_2m, O_1m}]
#   Reads the cached corrected-FS correlator matrix (fs_channels_v2), does the fixed-t0 GEVP with binsize
#   jackknife, and plots the effmass spectrum.  The corrected matrix is connected (tadpole diagrams vanish),
#   so the metric is no longer tadpole-degenerate.
#   Run:  ENS=... T0=3 BINSIZE=10 TSUB=0 python3 fs_channels_v2_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import glob
import numpy as np
import distill_contract_claude as dc

T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
TSUB = int(os.environ.get("TSUB", "0"))            # 1 = per-config t-sum (DC) subtraction before GEVP
SPLIT = int(os.environ.get("SPLIT", "1"))
NOP = int(os.environ.get("NOP", "3"))              # 2 -> drop O_1m ; 3 -> full
M2PS = 0.644


def gevp(Ct, C0):
    C0 = 0.5 * (C0 + C0.T)
    w, U = np.linalg.eigh(C0)
    keep = w > 1e-9 * w.max()
    Uk = U[:, keep] / np.sqrt(w[keep])
    M = Uk.T @ (0.5 * (Ct + Ct.T)) @ Uk
    return np.sort(np.linalg.eigvals(M).real)[::-1], int(keep.sum())


def main():
    tag = dc.ENS.split("nu0")[0]
    cands = glob.glob("fs_channels_v2_cache_claude/fs_channels_v2_%s_*cfg_d%d_claude.npy"
                      % (tag.replace(".", "p"), SPLIT))
    if not cands:
        print("# no cache found; run tmp_claude.sh first")
        return
    cache = max(cands, key=os.path.getsize)
    allC = np.load(cache)                            # (ncfg,3,3,twin)
    ncfg = allC.shape[0]
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    allC = allC[:, :NOP, :NOP, :]
    twin = allC.shape[-1]
    print("# GEVP on %s  ncfg=%d NOP=%d T0=%d BINSIZE=%d TSUB=%d" % (cache, ncfg, NOP, T0, BINSIZE, TSUB))
    if TSUB:
        allC = allC - allC[:, :, :, 1:].mean(axis=3, keepdims=True)

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    # metric conditioning report on the full mean
    Cm = blk.mean(0)
    print("# metric C(T0=%d) eigenvalues: %s" % (T0, np.array2string(np.linalg.eigvalsh(Cm[:, :, T0]), precision=4)))

    def effmass(Cmat):
        tmax = Cmat.shape[-1]
        ev = np.full((tmax, NOP), np.nan)
        for dt in range(tmax):
            if np.any(~np.isfinite(Cmat[:, :, dt])):
                continue
            try:
                e, nk = gevp(Cmat[:, :, dt], Cmat[:, :, T0])
                ev[dt, :len(e)] = e
            except Exception:
                pass
        with np.errstate(all="ignore"):
            return np.log(ev[:-1] / ev[1:])

    em_c = effmass(blk.mean(0))
    ems = np.array([effmass(np.delete(blk, i, 0).mean(0)) for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]

    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NOP)))
    for t in range(T0, min(tmax, 24)):
        row = "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_err[t, n]) for n in range(NOP))
        print("#  %2d | %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.7)
    ax.text(tmax * 0.6, M2PS + 0.012, r"$2m_{PS}\approx0.644$", fontsize=9, color="gray")
    cols = ["tab:green", "tab:red", "tab:blue"]
    mkr = ["o", "s", "^"]
    for n in range(NOP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.25)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n], marker=mkr[n], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.set_ylim(0.2, 1.1)
    ax.set_xlim(T0, min(tmax, 22))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("Corrected FS %d-op GEVP {$\\sigma^2_{00}$,$O_{2m}$,$O_{1m}$}  T0=%d  %s L1 %d cfg"
                 % (NOP, T0, tag, ncfg), fontsize=11)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_channels_v2_gevp_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
