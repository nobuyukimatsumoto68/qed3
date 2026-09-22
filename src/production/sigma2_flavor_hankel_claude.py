#!/usr/bin/env python3
# sigma2_flavor_hankel_claude.py  [chunk 1b: GEVP + Hankel+rebase on the {PP,FF,FP} flavor matrix at sigma^2_00]
#   The flavor ops PP=sigma_PS^2, FF=sigma_FS^2, FP=sigma_FS sigma_PS overlap the SAME 0++ tower with nearly
#   ORTHOGONAL amplitudes (<PP FF> tiny vs <PP PP>), so {PP,FF,FP} is a genuine variational basis (unlike the
#   redundant N_v-smearing copies).  Canonical Hankel+rebase (Dt=[0,2,4] reb NKEEP@REBT T0).  GEOM=0 (sigma^2_00).
#   Run: OFFSETS=0,2,4 REBT=4 NKEEP=2 T0=3 BINSIZE=10 python3 sigma2_flavor_hankel_claude.py

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
GEOM = int(os.environ.get("GEOM", "0"))
SPLIT = int(os.environ.get("SPLIT", "1"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,2,4").split(",")]
OPS = os.environ.get("OPS", "0,1,2")               # subset of {0=PP,1=FF,2=FP}
FLABS = ["PP", "FF", "FP"]
M2PS = 0.644


def hankel_reb(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    cache = max(glob.glob("sigma2_flavor_cache_claude/sigma2_flavor_%s_*cfg_g%d_d%d_claude.npy"
                          % (tag.replace(".", "p"), GEOM, SPLIT)), key=os.path.getsize)
    allC = np.load(cache)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    sel = [int(x) for x in OPS.split(",")]
    allC = allC[:, sel][:, :, sel]
    ncfg = allC.shape[0]
    print("# %s ncfg=%d ops=%s GEOM=%d OFFSETS=%s reb%d@%d T0=%d"
          % (cache.split("/")[-1], ncfg, [FLABS[s] for s in sel], GEOM, OFFSETS, NKEEP, REBT, T0))
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = hankel_reb(blk.mean(0), None)
    ems = np.array([hankel_reb(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]
    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)) + "  [ref 0.46 / 0.62]")
    for t in range(T0, min(tmax, 22)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_err[t, n]) for n in range(NKEEP))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ax.axhline(M2PS, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(tmax * 0.55, M2PS + 0.008, r"$2m_{PS}=0.644$", fontsize=9, color="gray")
    ax.axhline(0.46, color="tab:green", ls=":", lw=1, alpha=0.5)
    cols = ["tab:green", "tab:red", "tab:blue"]
    mkr = ["o", "s", "^"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_err[:, n]) & (em_err[:, n] < 0.25)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_err[g, n], color=cols[n % 3], marker=mkr[n % 3], ms=5, lw=1.1,
                    capsize=2.5, label="state %d" % n)
    ax.axvline(REBT, color="k", ls=":", lw=0.9, alpha=0.35)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(tmax, 20))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"flavor {%s} Hankel+rebase  GEOM=%d Dt=%s reb%d@%d  %s L1 %dcfg"
                 % (",".join(FLABS[s] for s in sel), GEOM, OFFSETS, NKEEP, REBT, tag, ncfg), fontsize=10)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_flavor_hankel_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
