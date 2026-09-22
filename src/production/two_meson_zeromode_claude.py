#!/usr/bin/env python3
# two_meson_zeromode_claude.py
# Per-config zero-mode (vacuum) subtraction of the PS two-meson correlator:
#   Asub(t) = A(t) - (1/W) sum_t A(t)   (removes the per-config t-constant = the vacuum <sigma sigma>).
# Plot <Asub(t)> jackknifed, LINEAR scale.  A(t) = normal-ordered PS.PS four-point (translation-averaged).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

CONTACT = float(os.environ.get("CONTACT", "0.5"))


def main():
    de.CONTACT = CONTACT
    print("# ENS=%s  ncfg=%d  CONTACT=%.2f  per-config zero-mode subtraction" % (dc.ENS.split("nu0")[0], len(dc.KS), CONTACT))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    FP = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        fp = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                acc += (dc.W10 * de.diags_pair(Phi, tau, s, s + dt)).sum()
            fp[dt] = 2.0 * (acc / ns).real
        FP.append(fp)
    FP = np.array(FP)                       # (ncfg, twin)
    ncfg, twin = FP.shape

    # per-config zero-mode subtraction
    Asub = FP - FP.mean(axis=1, keepdims=True)     # subtract window time-average per config

    # jackknife <Asub(t)>
    samp = np.array([np.delete(Asub, i, 0).mean(0) for i in range(ncfg)])
    m = samp.mean(0)
    err = np.sqrt((ncfg - 1) * np.mean((samp - m) ** 2, 0))
    print("\n  t    <Asub(t)>       err          |C/err|")
    for t in range(0, twin):
        print("  %2d  %+11.4e  %10.3e   %6.2f" % (t, m[t], err[t], abs(m[t]) / err[t]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tag = dc.ENS.split("nu0")[0]
    ts = np.arange(0, twin)
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.errorbar(ts, m, yerr=err, color="tab:red", marker="o", ms=4, lw=1, capsize=2)
    ax.axhline(0.0, color="k", lw=0.7)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle A^\mathrm{sub}(t)\rangle$  (PS$\cdot$PS, zero-mode subtracted)")
    ax.set_title("PS two-meson, per-config zero-mode subtraction (%s, %d cfg, contact=%.2f)" % (tag, ncfg, CONTACT))
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_zeromode_c%.1f_%s_claude.png" % (CONTACT, tag)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
