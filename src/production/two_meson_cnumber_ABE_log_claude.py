#!/usr/bin/env python3
# two_meson_cnumber_ABE_log_claude.py
# RAW (no t-sum subtraction) per-diagram correlators A(-S_S), B(-T_S), E(C_S^2) with the c-number diagrams,
# translation-averaged, config jackknife, LOG scale.  Plots |C(dt)| (sign printed in the table).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import two_meson_cnumber_claude as cn

IDX = [0, 1, 4]                                  # A, B, E
NAME = {0: "A(-S_S)", 1: "B(-T_S)", 4: "E(C_S^2)"}
STYLE = {0: ("tab:red", "o"), 1: ("tab:blue", "s"), 4: ("tab:green", "^")}


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  RAW A,B,E correlators (c-number diags, log)" % (tag, len(dc.KS)))
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
    allD = np.array(allD)
    ncfg, _, twin = allD.shape

    def jk(C):
        m = C.mean(0)
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])
        return m, np.sqrt((ncfg - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    means = {}
    errs = {}
    for i in IDX:
        means[i], errs[i] = jk(allD[:, i, :])
    print("\n  dt " + "".join("   %-16s" % NAME[i] for i in IDX))
    for dt in range(0, 14):
        row = "  %2d " % dt
        for i in IDX:
            row += "  %+10.3e(%.1e)" % (means[i][dt], errs[i][dt])
        print(row)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for i in IDX:
        col, mk = STYLE[i]
        m = means[i]
        e = errs[i]
        ax.errorbar(dts, np.abs(m[dts]), yerr=e[dts], color=col, marker=mk, ms=4, lw=0.9, capsize=1.5,
                    label=r"%s $|C|$" % NAME[i])
    ax.set_yscale("log")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$|C(dt)|$  (raw, no t-subtraction)")
    ax.set_title("Raw A,B,E per-diagram correlators, c-number diags (%s, %d cfg)" % (tag, ncfg))
    ax.legend(fontsize=10)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_cnumber_ABE_log_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
