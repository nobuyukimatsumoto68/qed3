#!/usr/bin/env python3
# ps_vs_fs_sigma2_free_claude.py
#   Free-limit sigma^2_00 diagonal effmass: PS (sigma_PS^2) vs FS (sigma_FS^2), side by side.
#   These are the FULL summed bilinears (not a leg-split); at m=0 the GW identity forces PS == FS.
#   Both must show the SAME spectrum, with the one-meson ground m_sigma = 0.378 as the lowest state
#   (genuine single-particle contamination of sigma^2 via the nonzero <sigma sigma^2> triangle; no
#   furnishing removes it -- only a variational GEVP does).
#   Source: the free 9-op flavor x geometry cache; op index = flavor*3+geom (flavor 0=PP=PS,1=FF=FS,2=FP).
#   PS sigma^2_00 = [0,0], FS sigma^2_00 = [3,3].
#   Run: python3 ps_vs_fs_sigma2_free_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np

CACHE = "sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_free_1cfg_d1_claude.npy"


def main():
    allC = np.load(CACHE)
    allC = 0.5 * (allC + np.swapaxes(allC, 1, 2))         # symmetrize op-ordering
    C = allC[0]                                            # single free config, (9,9,DT)
    DT = C.shape[-1]
    cPS = C[0, 0, :].real                                 # PP-s2 = sigma_PS^2 diagonal sigma^2_00
    cFS = C[3, 3, :].real                                 # FF-s2 = sigma_FS^2 diagonal sigma^2_00

    with np.errstate(all="ignore"):
        mPS = np.log(cPS[:-1] / cPS[1:])
        mFS = np.log(cFS[:-1] / cFS[1:])

    print("# sigma^2_00 diagonal effmass: PS vs FS (free L1; should be identical, both -> 0.378)")
    print("# refs: one-meson 0.378 ; (2,2) 0.556 ; two-meson 0.756")
    print("#  dt |  C_PS         m_PS   |  C_FS         m_FS   | |m_PS-m_FS|")
    for dt in range(1, DT - 1):
        mp = mPS[dt] if np.isfinite(mPS[dt]) else np.nan
        mf = mFS[dt] if np.isfinite(mFS[dt]) else np.nan
        d = abs(mp - mf) if (np.isfinite(mp) and np.isfinite(mf)) else np.nan
        print("#  %2d | %+.4e  %7.4f |  %+.4e  %7.4f |  %.2e" % (dt, cPS[dt], mp, cFS[dt], mf, d))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(DT - 1)
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for y, lab, col in [(0.378, r"one-meson $m_\sigma=0.378$", "tab:green"),
                        (0.556, r"$(2,2)=\{1,1,1,1\}$", "tab:orange"),
                        (0.756, r"two-meson $2m_\sigma=0.756$", "tab:blue")]:
        ax.axhline(y, color=col, ls="--", lw=1, alpha=0.55)
        ax.text(DT * 0.58, y + 0.008, lab, fontsize=9, color=col)
    gP = np.isfinite(mPS) & (cPS[:-1] * cPS[1:] > 0)
    gF = np.isfinite(mFS) & (cFS[:-1] * cFS[1:] > 0)
    ax.plot(ts[gP], mPS[gP], color="tab:red", marker="o", ms=6, lw=1.2, label=r"PS: $\sigma_{PS}^2$")
    ax.plot(ts[gF], mFS[gF], color="k", marker="x", ms=6, lw=0, label=r"FS: $\sigma_{FS}^2$ (overlaid)")
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(0, DT - 2)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE $\sigma^2_{00}$ effmass: PS vs FS (identical at $m=0$, GW) $\to$ one-meson $0.378$")
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/ps_vs_fs_sigma2_free_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
