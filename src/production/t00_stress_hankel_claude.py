#!/usr/bin/env python3
# t00_stress_hankel_claude.py  [CHUNK 1c: block-Hankel + rebase GEVP on the single T_00 correlator]
#
# The temporal-displacement basis {Dt=1,2,3} is COLLINEAR (rank-1 matrix -> degenerate GEVP), so it adds
# no variational power.  The team's standard way to sharpen a plateau from ONE operator's time series is
# the block-Hankel / GPOF GEVP: from the single folded correlator C(t) build Big(t)[a,b]=C(t+off[a]+off[b]),
# rebase, and take the moving-GEVP effmass.  We reuse the FROZEN core (do NOT fork its algebra):
#   final/analysis_axial/effmass_axial_tp_l3_perm_hankel_claude.py :: hankel_effmass_scalar.
#
# CANONICAL params (sigma/distillation thread, via qed3-7f 2026-09-16): offsets Dt=[0,2,4], reb2@4
# (NKEEP=2 @ REBT=4), T0=3, bin10.  Method refs: block-Hankel/GPOF Aubin-Orginos arXiv:1010.0202;
# GEVP Luscher-Wolff.  Units: em = log(lam_t/lam_{t+1}) dimensionless a_t*m (divide by a_t only at the end).
#
# Sign: the connected fermion-loop correlator here is NEGATIVE (overall -1); feed C_pos = -C (positive-
#   decaying) so the metric Big[t0] is positive-definite.  Free-limit = single complete-basis config, exact.
#
# Run:  ENS=free NVDIR=distill_Nv24 LREF=1 DTMAX=28 python3 t00_stress_hankel_claude.py
# See t00_stress_impl_plan_claude.md.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "final/analysis_axial")
import numpy as np
import distill_contract_claude as dc
import t00_stress_claude as t0s
import effmass_axial_tp_l3_perm_hankel_claude as hk

DTMAX = int(os.environ.get("DTMAX", "28"))
OFFS = [int(x) for x in os.environ.get("OFFSETS", "0,2,4").split(",")]
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "2"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "1"))
NCFG = int(os.environ.get("NCFG", "0"))
AT = float(os.environ.get("AT", "0.2"))
ROVER = 1.0 / 0.189


def corr_all(ks, w):
    # single-operator T_00 connected correlator C(dt) per config (dt=0..DTMAX-1), sign-flipped positive.
    t0s.DTMAX = DTMAX
    out = []
    for k in ks:
        C, _ = t0s.corr_one_config(k, w)
        out.append(-C)                              # positive-decaying
    return np.array(out)                            # (ncfg, DTMAX)


def main():
    tag = dc.ENS
    dual = dc.dual_areas_from_mesh()
    w = dual * dc.Y00
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    if not ks:
        print("# no perambulators found in %s" % dc.PERAM_DIR)
        return
    print("# ENS=%s NVDIR=%s  offsets=%s reb%d@%d T0=%d bin%d  ncfg=%d"
          % (tag, dc.NVDIR, OFFS, NKEEP, REBT, T0, BINSIZE, len(ks)))

    allC = corr_all(ks, w)                           # (ncfg, DTMAX), already sign-flipped positive
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])

    em_c, Vop = hk.hankel_effmass_scalar(blk.mean(0), OFFS, REBT, NKEEP, T0, AT)
    if nb < 2:
        em_err = np.zeros_like(em_c)                 # single config (free limit): exact, no jackknife
    else:
        ems = []
        for i in range(nb):
            Ci = np.delete(blk, i, 0).mean(0)
            em_i, _ = hk.hankel_effmass_scalar(Ci, OFFS, REBT, NKEEP, T0, AT, Vop=Vop)
            ems.append(em_i)
        ems = np.array(ems)
        em_err = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, axis=0))

    tmax1 = em_c.shape[0]
    print("\n#  t |  ground(err)     1st-exc(err)   [free refs sigma=%.3f T00=%.3f 2m=%.3f]"
          % (2.0 / ROVER, 3.0 / ROVER, 4.0 / ROVER))
    for t in range(tmax1):
        g0 = em_c[t, 0]
        e0 = em_err[t, 0]
        g1 = em_c[t, 1] if NKEEP > 1 else np.nan
        e1 = em_err[t, 1] if NKEEP > 1 else np.nan
        if not np.isfinite(g0):
            continue
        print("#  %2d | %7.4f(%.4f)   %7.4f(%.4f)" % (t, g0, e0, g1, e1))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tt = np.arange(tmax1)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R$", "tab:green"),
                          (3.0 / ROVER, r"$T_{00}=3/R$", "tab:red"),
                          (4.0 / ROVER, r"$2m=4/R$", "gray")):
        ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
        ax.text(0.2, val + 0.006, lab, fontsize=9, color=col)
    g = np.isfinite(em_c[:, 0]) & (em_err[:, 0] < 0.3)
    ax.errorbar(tt[g], em_c[g, 0], yerr=em_err[g, 0], color="tab:red", marker="o", ms=5, lw=1.1,
                capsize=2.5, label="ground (Hankel)")
    if NKEEP > 1:
        ge = np.isfinite(em_c[:, 1]) & (em_err[:, 1] < 0.3)
        ax.errorbar(tt[ge], em_c[ge, 1], yerr=em_err[ge, 1], color="tab:blue", marker="s", ms=5, lw=1.0,
                    capsize=2.5, markerfacecolor="none", ls=":", label="1st-exc (Hankel)")
    ax.set_ylim(0.2, 1.1)
    ax.set_xlim(0, min(tmax1, 20))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}=\log(\lambda_t/\lambda_{t+1})$")
    ax.set_title(r"$T_{00}$ block-Hankel  off%s reb%d@%d T0=%d  %s  %d cfg"
                 % ("-".join(map(str, OFFS)), NKEEP, REBT, T0, tag, ncfg), fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_stress_hankel_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
