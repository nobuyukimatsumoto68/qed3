#!/usr/bin/env python3
# sigma2_ppff_cross_effmass_claude.py
#   Effective mass of the PP-FF cross channel of the sigma^2 flavor-geometry correlator.
#   Hypothesis (NM): the two-meson P^2 and F^2 states are DIFFERENT states that mix only through the
#   single-meson (2,2); hence the PP-FF cross block carries ONLY the single-meson pole (not the
#   two-meson), while the PP/FF diagonals carry both.  This driver isolates the cross channel two ways:
#     (a) leading singular value s0(dt) of the 3x3 cross block  B(dt) = C[PP-geom, FF-geom](dt)
#     (b) the sigma^2_00-sigma^2_00 element  C[0,3](dt)
#   and effmasses each (plateau-subtracted, single free config).  Compare to the single-meson (2,2)
#   ref and the two-meson 2m_PS ref.
#   Run: LTAG=L1 python3 sigma2_ppff_cross_effmass_claude.py
#        LTAG=L2 CACHE=... M_PS=0.393 STATE_A=0.69 TWO_MESON=0.786 python3 sigma2_ppff_cross_effmass_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np

LTAG = os.environ.get("LTAG", "L1")
CACHE = os.environ.get("CACHE",
    "sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_free_1cfg_d1_claude.npy")
STATE_1111 = float(os.environ.get("STATE_A", "0.556"))   # single-meson (2,2): L1 0.556, L2 ~0.69
TWO_MESON = float(os.environ.get("TWO_MESON", "0.756"))  # 2 m_PS: L1 0.756, L2 0.786
PLAT_LO = int(os.environ.get("PLAT_LO", "16"))           # plateau region for the disconnected subtraction
DTMAX = int(os.environ.get("DTMAX", "16"))


def effmass_plat(c):
    # connected: subtract the large-dt plateau (disconnected <PP><FF>), then log ratio
    plat = c[PLAT_LO:].mean()
    cc = c - plat
    with np.errstate(all="ignore"):
        return np.log(cc[:-1] / cc[1:])


def main():
    C = np.load(CACHE)
    C = 0.5 * (C + np.swapaxes(C, 1, 2))
    M = C.mean(0)                                        # (9,9,twin)
    twin = M.shape[-1]

    # cross block B(dt) = C[PP-geom (0,1,2), FF-geom (3,4,5)]
    s0 = np.zeros(twin)
    for dt in range(twin):
        B = M[0:3, 3:6, dt]
        s0[dt] = np.linalg.svd(B, compute_uv=False)[0]  # leading singular value
    c_s2 = M[0, 3, :]                                    # sigma^2_00 PP-FF element

    em_s0 = effmass_plat(s0)
    em_s2 = effmass_plat(c_s2)
    # diagonal PP sigma^2_00 for contrast (carries both poles)
    em_diag = effmass_plat(M[0, 0, :])

    print("# FREE %s  PP-FF CROSS channel effmass  (single free config)" % LTAG)
    print("# refs: single-meson (2,2) = %.3f ;  two-meson 2m_PS = %.3f" % (STATE_1111, TWO_MESON))
    print("#  dt |  cross SVD0   cross s2-elt |  PP diag s2 (both poles)")
    for dt in range(1, DTMAX):
        def f(x):
            return "% .4f" % x if np.isfinite(x) else "  ---  "
        print("#  %2d |  %s     %s   |   %s" % (dt, f(em_s0[dt]), f(em_s2[dt]), f(em_diag[dt])))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, DTMAX)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.axhline(STATE_1111, color="tab:orange", ls="--", lw=1, alpha=0.7)
    ax.text(DTMAX * 0.5, STATE_1111 + 0.008, r"single-meson (2,2) $=%.3f$" % STATE_1111,
            color="tab:orange", fontsize=9)
    ax.axhline(TWO_MESON, color="tab:blue", ls="--", lw=1, alpha=0.7)
    ax.text(DTMAX * 0.5, TWO_MESON + 0.008, r"two-meson $2m_{PS}=%.3f$" % TWO_MESON,
            color="tab:blue", fontsize=9)
    # cross channel: red circle ; diagonal (both poles): black square (color-blind safe)
    ax.plot(dts, em_s0[dts], color="tab:red", marker="o", ms=6, lw=1.2, label="PP-FF cross (leading SVD)")
    ax.plot(dts, em_s2[dts], color="tab:purple", marker="^", ms=5, lw=1.0, ls="--",
            label=r"PP-FF cross ($\sigma^2_{00}$ elt)")
    ax.plot(dts, em_diag[dts], color="black", marker="s", ms=5, lw=1.0, alpha=0.7,
            label=r"PP diagonal $\sigma^2_{00}$ (both poles)")
    ax.set_ylim(0.3, 1.2)
    ax.set_xlim(1, DTMAX)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE %s  PP-FF cross-channel effmass  (single-meson only?)" % LTAG)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_ppff_cross_effmass_%s_claude.png" % LTAG
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
