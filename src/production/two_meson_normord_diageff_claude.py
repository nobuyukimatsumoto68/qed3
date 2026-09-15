#!/usr/bin/env python3
# two_meson_normord_diageff_claude.py
# Per-diagram RAW effmass (log C(dt)/C(dt+1), no plateau subtraction) for ALL A--J with the
# COLLISION-contact-removed diagrams (:sigma^2: = diag_effmass.diags_pair, tau(s,s)-1/2 I everywhere).
# Config jackknife.  Shows which diagrams plateau above 0.6 (two-sigma/0++) and which do not.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

LABELS = de.LABELS
at = 0.2
atm_f2 = at * 3.08


def main():
    tag = dc.ENS.split("nu0")[0]
    de.CONTACT = 0.5                                  # collision contact removed everywhere
    print("# ENS=%s  ncfg=%d  per-diagram RAW effmass, :sigma^2: (tau-1/2 I everywhere)" % (tag, len(dc.KS)))
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
    allD = np.array(allD)
    ncfg, _, twin = allD.shape
    allD = allD - allD.mean(axis=2, keepdims=True)   # per-config t-total (time-average) subtraction

    def jk_corr(C):                                  # jk mean, err of the (t-sub) correlator
        m = C.mean(0)
        samp = np.array([np.delete(C, i, 0).mean(0) for i in range(ncfg)])
        return m, np.sqrt((ncfg - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))

    means = {}
    print("\n  diagram        C(dt2)      C(dt5)      C(dt8)   (t-total subtracted)")
    for i in range(10):
        m, e = jk_corr(allD[:, i, :])
        means[i] = (m, e)
        print("  %-14s  %+9.2e  %+9.2e  %+9.2e" % (LABELS[i], m[2], m[5], m[8]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, 16)
    fig, axs = plt.subplots(2, 5, figsize=(20, 7.5))
    for i in range(10):
        ax = axs[i // 5, i % 5]
        m, e = means[i]
        ax.errorbar(dts, np.abs(m[dts]), yerr=e[dts], marker="o", ms=3, lw=0.8, capsize=1.5, color="tab:red")
        ax.set_yscale("log")
        ax.set_title(LABELS[i], fontsize=10)
        ax.set_xlabel("dt")
    fig.suptitle(":$\\sigma^2$: per-diagram |C(dt)|, collision-contact removed, t-total subtracted "
                 "(%s, %d cfg)  [log, jackknife]" % (tag, ncfg))
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_normord_diagcorr_%s_claude.png" % tag
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
