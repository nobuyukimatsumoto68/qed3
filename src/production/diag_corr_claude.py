#!/usr/bin/env python3
# diag_corr_claude.py -- just the per-diagram two-meson correlators <A_i(dt)>, LINEAR, jackknifed.
# Translation-averaged over window sources. NO plateau subtraction. CONTACT env (0=raw, 0.5=normal-ordered).
# Shows which diagrams have a vacuum (flat, nonzero at large dt) vs decay to zero.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

CONTACT = float(os.environ.get("CONTACT", "0.0"))
LABELS = de.LABELS


def main():
    de.CONTACT = CONTACT
    print("# ENS=%s  ncfg=%d  CONTACT=%.2f  per-diagram correlators" % (dc.ENS.split("nu0")[0], len(dc.KS), CONTACT))
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
                acc += de.diags_pair(Phi, tau, s, s + dt)
            D[:, dt] = (acc / ns).real
        allD.append(D)
    allD = np.array(allD)                 # (ncfg, 10, twin)
    ncfg, _, twin = allD.shape
    # subtract the per-config time-average over the window (the vector-disc DC/zero-mode subtraction)
    allD = allD - allD.mean(axis=2, keepdims=True)

    def jk(C):                            # C (ncfg, twin) -> mean, err
        m = C.mean(0)
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])
        return m, np.sqrt((ncfg - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tag = dc.ENS.split("nu0")[0]
    dts = np.arange(1, 16)
    fig, axs = plt.subplots(2, 5, figsize=(20, 7.5))
    for i in range(10):
        ax = axs[i // 5, i % 5]
        m, e = jk(allD[:, i, :])
        ax.errorbar(dts, m[dts], yerr=e[dts], marker="o", ms=3, lw=0.8, capsize=1.5, color="tab:red")
        ax.axhline(0.0, color="k", lw=0.6)
        ax.set_title(LABELS[i], fontsize=10)
        ax.set_xlabel("dt")
    fig.suptitle("Per-diagram two-meson correlators, per-config time-avg subtracted <A_i^sub(dt)> "
                 "(%s, %d cfg, contact=%.2f)  [linear, jackknife]" % (tag, ncfg, CONTACT))
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs("figs", exist_ok=True)
    out = "figs/diag_corr_dcsub_c%.1f_%s_claude.png" % (CONTACT, tag)
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("# -> %s" % out)


if __name__ == "__main__":
    main()
