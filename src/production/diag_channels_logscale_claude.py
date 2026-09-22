#!/usr/bin/env python3
# diag_channels_logscale_claude.py
#   Diagnostic: the 3 diagonal geometry correlators {sigma^2_00, O_2m, O_1m} of a free cache, |C| in LOG
#   scale, to see where the single-config signal hits the solve-tolerance noise floor (why the GEVP is junk).
#   Sign flips are marked (open markers = negative C).
#   Run: CACHE=... LTAG=L2 python3 diag_channels_logscale_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import numpy as np

LTAG = os.environ.get("LTAG", "L2")
CACHE = os.environ.get("CACHE", "sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_free_L2_1cfg_d1_claude.npy")
LABELS = ["sigma^2_00", "O_2m", "O_1m"]


def main():
    C = np.load(CACHE)[0]                       # (9,9,DT)
    DT = C.shape[-1]
    diag = [C[i, i, :].real for i in range(3)]  # PP geometry diagonals

    print("# %s diagonal channels (raw, incl. sign)" % LTAG)
    print("#  t |   sigma^2_00        O_2m             O_1m")
    for t in range(1, DT):
        print("#  %2d | %+.5e   %+.5e   %+.5e" % (t, diag[0][t], diag[1][t], diag[2][t]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(DT)
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    cols = ["tab:red", "tab:blue", "tab:green"]
    mkr = ["o", "s", "^"]
    for i in range(3):
        y = diag[i]
        pos = y > 0
        neg = y < 0
        ax.plot(ts[pos], np.abs(y[pos]), color=cols[i], marker=mkr[i], ms=5, lw=1.0, label=LABELS[i] + " (+)")
        ax.plot(ts[neg], np.abs(y[neg]), color=cols[i], marker=mkr[i], ms=7, lw=0, mfc="none", label=LABELS[i] + " (-)")
    ax.set_yscale("log")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|C(t)|$")
    ax.set_title(r"FREE %s diagonal channels $\{\sigma^2_{00},O_{2m},O_{1m}\}$ (log; open = negative)" % LTAG)
    ax.legend(fontsize=8, ncol=3)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_channels_logscale_%s_claude.png" % LTAG
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
