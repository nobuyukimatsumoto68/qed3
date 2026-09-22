#!/usr/bin/env python3
# fs_gevp_3x3_free_plot_claude.py
#   FREE-limit 3x3 geometry GEVP {sigma^2_00, O_2m (antipodal), O_1m (coincident split)} effmass plot.
#   Reads the cached free matrix built by fs_gevp_point_claude.run_gevp (single exact config).  The plain
#   run_gevp goes nan because C0 at T0 is indefinite (O_1m has a negative small-dt piece) -> the rank-drop
#   makes a 3-vs-2 shape mismatch.  Here we take the FULL generalized eigenvalues eig(C0^{-1} C(t)) (no
#   state dropped) so all three levels are kept, then effmass = log(lam(t)/lam(t+1)).
#   Reference free scales: one-meson m_sigma = 2E0 = 0.378 ; (2,2) single meson 2E1 ~ 0.556 ; two-meson
#   2 m_sigma = 0.756.
#   Run:  python3 fs_gevp_3x3_free_plot_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np

T0 = int(os.environ.get("T0", "3"))
CACHE = "fs_gevp_cache_claude/fs_gevp_point_free_1cfg_d1_claude.npy"


def gevp_full(Ct, C0):
    # full generalized eigenvalues, no state dropped ; C0 may be indefinite
    M = np.linalg.solve(C0, Ct)
    ev = np.linalg.eigvals(M)
    return np.sort(ev.real)[::-1]


def main():
    C = np.load(CACHE)[0]
    C = 0.5 * (C + np.swapaxes(C, 0, 1))
    DT = C.shape[-1]
    C0 = 0.5 * (C[:, :, T0] + C[:, :, T0].T)

    lam = np.full((DT, 3), np.nan)
    for dt in range(DT):
        if np.any(~np.isfinite(C[:, :, dt])):
            continue
        Ct = 0.5 * (C[:, :, dt] + C[:, :, dt].T)
        lam[dt] = gevp_full(Ct, C0)

    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])

    print("#  t |   m0        m1        m2      [one-meson 0.378 / (2,2) 0.556 / two-meson 0.756]")
    for t in range(T0 + 1, min(DT - 1, 20)):
        row = "  ".join("%7.4f" % em[t, i] if np.isfinite(em[t, i]) else "  ---  " for i in range(3))
        print("#  %2d | %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(em.shape[0])
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for y, lab, col in [(0.378, r"one-meson $m_\sigma=0.378$", "tab:green"),
                        (0.556, r"$(2,2)=\{1,1,1,1\}$ $2E_1$", "tab:orange"),
                        (0.756, r"two-meson $2m_\sigma=0.756$", "tab:blue")]:
        ax.axhline(y, color=col, ls="--", lw=1, alpha=0.6)
        ax.text(DT * 0.60, y + 0.008, lab, fontsize=9, color=col)
    cols = ["tab:green", "tab:orange", "tab:blue"]
    mkr = ["o", "s", "^"]
    for n in range(3):
        g = np.isfinite(em[:, n]) & (lam[:-1, n] * lam[1:, n] > 0)
        ax.plot(ts[g], em[g, n], color=cols[n], marker=mkr[n], ms=5, lw=1.1, label="state %d" % n)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(T0, min(DT - 1, 20))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE 3x3 geometry GEVP $\{\sigma^2_{00}, O_{2m}, O_{1m}\}$  (T0=%d)" % T0)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/fs_gevp_3x3_free_effmass_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
