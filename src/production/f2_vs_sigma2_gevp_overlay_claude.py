#!/usr/bin/env python3
# f2_vs_sigma2_gevp_overlay_claude.py
#   Overlay the purely-GLUONIC F^2 effmass (glueball GEVP) on the purely-FERMIONIC {sigma^2_00, O_2m, O_1m}
#   3-op GEVP spectrum, in COMMON physical units (m = a_t m / a_t).  NO coupled matrix (that 2x2 was
#   ill-conditioned) -- just the two independent spectra side by side to see where F^2 (2.88) sits relative
#   to the sigma^2 one-meson-rich ground / two-meson.
#     - fermionic: fs_gevp_point {sigma^2_00,O_2m,O_1m} 3x3 GEVP (cache), lattice effmass log(ev[t]/ev[t+1])
#       times 1/a_t  -> physical.  x = timeslice t.
#     - gluonic: glueball GEVP Delta_eff = -log(lambda)/(dt*a_t) (already physical), state 0 = F^2.
#       x = t_phys/a_t = timeslice.
#   Run: python3 f2_vs_sigma2_gevp_overlay_claude.py [glue_dat]

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

AT = float(os.environ.get("AT", "0.2"))
T0 = int(os.environ.get("T0", "3"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
DTMAX = int(os.environ.get("DTMAX", "24"))
GLUE_DAT = sys.argv[1] if len(sys.argv) > 1 else \
    "/tmp/claude-1000/-mnt-barracuda22-qed3/cdb76ded-4ed0-4f78-8eda-57542aa97f37/scratchpad/f2_gevp_repro.dat"


def gevp(Ct, C0):
    C0 = 0.5 * (C0 + C0.T)
    w, Uv = np.linalg.eigh(C0)
    keep = w > 1e-10 * w.max()
    Uk = Uv[:, keep] / np.sqrt(w[keep])
    M = Uk.T @ (0.5 * (Ct + Ct.T)) @ Uk
    return np.sort(np.linalg.eigvals(M).real)[::-1]


def effmass(Cm, nstate):
    ev = np.full((DTMAX, nstate), np.nan)
    for dt in range(DTMAX):
        if not np.all(np.isfinite(Cm[:, :, dt])):
            continue
        try:
            e = gevp(Cm[:, :, dt], Cm[:, :, T0])
            m = min(nstate, len(e))
            ev[dt, :m] = e[:m]                    # metric may prune to rank<nstate
        except Exception:
            pass
    with np.errstate(all="ignore"):
        r = ev[:-1] / ev[1:]
        r[r <= 0] = np.nan                       # drop sign-flips (state went negative)
        return np.log(r)                         # lattice a_t m


def main():
    tag = dc.ENS.split("nu0")[0]
    # --- fermionic 3-op GEVP from cache ---
    cache = "fs_gevp_cache_claude/fs_gevp_point_%s_400cfg_d1_claude.npy" % tag.replace(".", "p")
    allC = np.load(cache)                          # (ncfg,3,3,DTMAX)
    ncfg = allC.shape[0]
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    ns = 2                                         # 3-op metric prunes to rank 2 at T0
    em_c = effmass(blk.mean(0), ns)
    ems = np.array([effmass(np.delete(blk, i, 0).mean(0), ns) for i in range(nb)])
    em_err = np.sqrt((nb - 1) * np.nanmean((ems - np.nanmean(ems, 0)) ** 2, axis=0))
    em_phys = em_c / AT                            # -> physical
    em_phys_err = em_err / AT
    print("# fermionic {sigma^2_00,O_2m,O_1m} GEVP (physical m = a_t m / a_t; rank-2):")
    print("#  t | state0(err)     state1(err)")
    for t in range(T0, min(DTMAX - 1, 18)):
        print("#  %2d | %6.3f(%.3f)  %6.3f(%.3f)"
              % (t, em_phys[t, 0], em_phys_err[t, 0], em_phys[t, 1], em_phys_err[t, 1]))

    # --- gluonic F^2 from glueball GEVP .dat (col: t ground_m ground_e s0_m s0_e s1_m s1_e ; s0=F^2) ---
    g = np.loadtxt(GLUE_DAT)
    gt_slice = np.rint(g[:, 0] / AT).astype(int)   # physical t -> timeslice
    gF = g[:, 3]                                    # state 0 = F^2 (already physical)
    gFe = g[:, 4]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(DTMAX - 1)
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    cols = ["tab:blue", "tab:red", "tab:green"]
    mks = ["o", "s", "^"]
    labs = [r"$\sigma^2$ GEVP state 0 (1-meson-rich)", r"state 1 (two-meson)"]
    for n in range(2):
        m = np.isfinite(em_phys[:, n]) & np.isfinite(em_phys_err[:, n]) & (em_phys_err[:, n] < 1.5)
        ax.errorbar(ts[m], em_phys[m, n], yerr=em_phys_err[m, n], color=cols[n], marker=mks[n],
                    ms=5, lw=1.1, capsize=2.5, label=labs[n])
    gm = gFe < 1.0
    ax.errorbar(gt_slice[gm], gF[gm], yerr=gFe[gm], color="black", marker="D", ms=6, lw=1.4,
                capsize=3, label=r"$F^2$ glueball GEVP (state 0)")
    ax.axhline(2.88, color="gray", ls="--", lw=0.9, alpha=0.6)
    ax.text(13.0, 2.92, r"$F^2$ fit $2.88$", fontsize=8, color="gray")
    ax.set_ylim(0.8, 3.6)
    ax.set_xlim(0.5, 17)
    ax.set_xlabel(r"timeslice $t$")
    ax.set_ylabel(r"physical mass  $m = a_t m_{\rm eff}/a_t$")
    ax.set_title(r"$F^2$ vs $\{\sigma^2_{00},O_{2m},O_{1m}\}$ GEVP  %s L1 %d cfg (bin %d, $T_0$=%d)"
                 % (tag, ncfg, BINSIZE, T0), fontsize=11)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/f2_vs_sigma2_gevp_overlay_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
