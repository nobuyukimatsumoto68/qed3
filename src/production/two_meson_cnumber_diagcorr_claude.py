#!/usr/bin/env python3
# two_meson_cnumber_diagcorr_claude.py
# Per-diagram two-meson correlators with the CORRECT c-number (tadpole-only) diagrams
# (`two_meson_cnumber_claude.diags_pair_cnumber`): raw tau in the connected pieces, D_S -> D_S - c0.
# Per-config t-sum (time-average over the window) subtracted, then config jackknife; LINEAR scale.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import two_meson_cnumber_claude as cn
import diag_effmass_claude as de       # LABELS


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  per-diagram corr (c-number diags, t-sum subtracted, linear)"
          % (tag, len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    allD = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
        D = np.zeros((10, twin))
        for dt in range(twin):
            ns = twin - dt
            acc = np.zeros(10, complex)
            for s in range(ns):
                acc += cn.diags_pair_cnumber(Phi, tau, s, s + dt)
            D[:, dt] = (acc / ns).real
        allD.append(D)
    allD = np.array(allD)                       # (ncfg, 10, twin)
    ncfg, _, twin = allD.shape
    # per-config t-sum (time-average) subtraction
    allD = allD - allD.mean(axis=2, keepdims=True)

    def jk(C):                                  # C (ncfg, twin) -> mean, err
        m = C.mean(0)
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])
        return m, np.sqrt((ncfg - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, 16)
    fig, axs = plt.subplots(2, 5, figsize=(20, 7.5))
    for i in range(10):
        ax = axs[i // 5, i % 5]
        m, e = jk(allD[:, i, :])
        ax.errorbar(dts, m[dts], yerr=e[dts], marker="o", ms=3, lw=0.8, capsize=1.5, color="tab:red")
        ax.axhline(0.0, color="k", lw=0.6)
        ax.set_title(de.LABELS[i], fontsize=10)
        ax.set_xlabel("dt")
    fig.suptitle("Per-diagram two-meson correlators, c-number diags, t-sum subtracted "
                 "(%s, %d cfg)  [linear, jackknife]" % (tag, ncfg))
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_cnumber_diagcorr_%s_claude.png" % tag
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("# -> %s" % out)
    # print a compact table of the subtracted correlator at a few dt
    print("\n  diagram        C_sub(dt=2)     C_sub(dt=5)     C_sub(dt=8)")
    for i in range(10):
        m, e = jk(allD[:, i, :])
        print("  %-14s  %+10.3e     %+10.3e     %+10.3e" % (de.LABELS[i], m[2], m[5], m[8]))


if __name__ == "__main__":
    main()
