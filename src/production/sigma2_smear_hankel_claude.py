#!/usr/bin/env python3
# sigma2_smear_hankel_claude.py  [chunk 2: GEVP + Hankel+rebase on the N_v-smearing-expanded PS^2 point basis]
#   Reads the 12-op (3 point ops x 4 smearings {4,8,16,24}) connected matrix from sigma2_smear_gevp_v2_claude,
#   symmetrizes at the ENSEMBLE level (O_1m is time-split -> raw matrix asymmetric), then runs the canonical
#   Hankel+rebase GEVP (same machinery as fs_channels_v2_hankel: Dt=OFFSETS, rebase NKEEP@REBT, metric T0).
#   Goal: does the multi-smearing variational basis sharpen the two-meson m1 / resolve a 3rd state vs the
#   3-op (Nv=24-only) result (state0~0.46, state1~0.62)?  Refs: GPOF Aubin-Orginos 1010.0202; Blossier 0902.1265.
#
#   OPS selects a column subset of the 12 (default all).  Run:
#     OFFSETS=0,2,4 REBT=4 NKEEP=2 T0=3 BINSIZE=10 python3 sigma2_smear_hankel_claude.py

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
SPLIT = int(os.environ.get("SPLIT", "1"))
SMEARS = [int(x) for x in os.environ.get("SMEARS", "4,8,16,24").split(",")]
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,2,4").split(",")]
NSM = len(SMEARS)
OPLAB = ["s2_00", "O_2m", "O_1m"]
# default: all 12 ops ; OPS env = comma list of flat indices (op*NSM + ismear)
OPS = os.environ.get("OPS", "")
M2PS = 0.644
AT = 0.2


def hankel_reb(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    cands = glob.glob("sigma2_smear_cache_claude/sigma2_smear_%s_*cfg_sm%s_d%d_claude.npy"
                      % (tag.replace(".", "p"), "-".join(map(str, SMEARS)), SPLIT))
    if not cands:
        print("# NO cache yet (matrix still building?) -- expected sigma2_smear_cache_claude/sigma2_smear_%s_*_sm%s_d%d"
              % (tag.replace(".", "p"), "-".join(map(str, SMEARS)), SPLIT))
        return
    cache = max(cands, key=os.path.getsize)
    allC = np.load(cache)                               # (ncfg, 12, 12, DTMAX)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))       # ENSEMBLE-level symmetrize (O_1m time-split)
    ncfg = allC.shape[0]
    if OPS:
        sel = [int(x) for x in OPS.split(",")]
        allC = allC[:, sel][:, :, sel]
    nop = allC.shape[1]
    labs = ["%s@%d" % (OPLAB[i // NSM], SMEARS[i % NSM]) for i in range(allC.shape[1])] if not OPS \
        else ["%s@%d" % (OPLAB[int(x) // NSM], SMEARS[int(x) % NSM]) for x in OPS.split(",")]
    print("# cache <- %s  ncfg=%d nop=%d  OFFSETS=%s REBT=%d NKEEP=%d T0=%d" % (cache, ncfg, nop, OFFSETS, REBT, NKEEP, T0))
    print("# ops: %s" % ", ".join(labs))

    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = hankel_reb(blk.mean(0), None)
    ems = np.array([hankel_reb(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]

    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)) + "  [ref: state0~0.46, state1~0.62]")
    for t in range(T0, min(tmax, 22)):
        row = "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_err[t, n]) for n in range(NKEEP))
        print("#  %2d | %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(tmax * 0.6, M2PS + 0.008, r"$2m_{PS}\approx0.644$", fontsize=9, color="gray")
    cols = ["tab:green", "tab:red", "tab:blue", "tab:purple"]
    mkr = ["o", "s", "^", "D"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.25)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 4], marker=mkr[n % 4], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(tmax, 20))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"PS $\sigma^2$ $N_v$-smear-expanded Hankel+rebase  nop=%d Dt=%s reb%d@%d T0=%d  %s L1 %dcfg"
                 % (nop, OFFSETS, NKEEP, REBT, T0, tag, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_smear_hankel_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
