#!/usr/bin/env python3
# fs_channels_v2_hankel_claude.py  [Hankel+rebase on the corrected FS 3x3 channel matrix]
#   Block-Hankel (offsets from OFFSETS) on the corrected-FS {sigma^2_00, O_2m, O_1m} matrix, then rebase
#   NKEEP states at REBT (metric T0).  Reuses hankel_off/staged_project/rebased_effmass_fixed.
#   Run:  OFFSETS=0,2,4 REBT=4 NKEEP=2 T0=3 BINSIZE=10 NOP=3 python3 fs_channels_v2_hankel_claude.py

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
import hankel_rebase_scan_claude as hs

T0 = int(os.environ.get("T0", "3"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "2"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NOP = int(os.environ.get("NOP", "3"))
SPLIT = int(os.environ.get("SPLIT", "1"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,2,4").split(",")]
M2PS = 0.644


def hankel_reb(Cmat, Vfix):
    # Cmat (NOP,NOP,twin) -> Hankel(OFFSETS) -> rebase (Vfix or build) -> effmass
    Cts = np.transpose(Cmat, (2, 0, 1))              # (twin, NOP, NOP)
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    cands = glob.glob("fs_channels_v2_cache_claude/fs_channels_v2_%s_*cfg_d%d_claude.npy"
                      % (tag.replace(".", "p"), SPLIT))
    cache = max(cands, key=os.path.getsize)
    allC = np.load(cache)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    allC = allC[:, :NOP, :NOP, :]
    ncfg = allC.shape[0]
    print("# Hankel+rebase FS  ncfg=%d NOP=%d OFFSETS=%s REBT=%d NKEEP=%d T0=%d BIN=%d"
          % (ncfg, NOP, OFFSETS, REBT, NKEEP, T0, BINSIZE))

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = hankel_reb(blk.mean(0), None)
    ems = np.array([hankel_reb(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]

    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)))
    for t in range(T0, min(tmax, 24)):
        row = "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_err[t, n]) for n in range(NKEEP))
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
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.25)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n], marker=mkr[n], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.4)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(tmax, 22))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("Corrected FS Hankel+rebase  Dt=%s NOP=%d rebase %d@%d T0=%d  %s L1 %d cfg"
                 % (OFFSETS, NOP, NKEEP, REBT, T0, tag, ncfg), fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_channels_v2_hankel_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
