#!/usr/bin/env python3
# diagramA_corr_linear_claude.py
#   Diagram-A sigma^2_00 correlator C_A(t) in LINEAR scale (S-leg, contact-subtracted), single free config.
#   Purpose: understand why diagram A's L2 effmass goes junk -- does the correlator dip negative / cross zero?
#   Run: ENS=free LREF=2 NVDIR=distill_Nv84 LTAG=L2 python3 diagramA_corr_linear_claude.py

import os
os.environ.setdefault("ENS", "free")
os.environ.setdefault("LREF", "2")
os.environ.setdefault("NVDIR", "distill_Nv84")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import fs_gevp_point_perclass_claude as pc

DTMAX = int(os.environ.get("DTMAX", "24"))
LTAG = os.environ.get("LTAG", "L2")


def main():
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(dc.KS[0])
    dual = dc.dual_areas_from_mesh().astype(float)
    wY = dual * dc.Y00
    Pmap = G.antipodal_map()
    cls_of = [pc.classify(cy) for cy in G.PERMS]

    CA = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        vt = [dt, dt, 0, 0]
        offs = set((vt[a], vt[b]) for a in range(4) for b in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        vspec = G.op_vspec(0, ('i', 'j'), dual, wY) + G.op_vspec(0, ('k', 'l'), dual, wY)
        acc = 0.0
        for ip, cyc in enumerate(G.PERMS):
            if cls_of[ip] == "A":
                acc += G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap).real
        CA[dt] = acc / len(s0s)

    print("# %s diagram-A sigma^2_00 correlator C_A(t) (linear)" % LTAG)
    print("#  dt |   C_A(dt)")
    for dt in range(1, DTMAX):
        if np.isfinite(CA[dt]):
            print("#  %2d | %+.8e" % (dt, CA[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(DTMAX)
    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    LOG = int(os.environ.get("LOG", "0"))
    if LOG:
        pos = CA > 0
        neg = CA < 0
        ax.plot(ts[pos], np.abs(CA[pos]), color="tab:red", marker="o", ms=5, lw=1.2, label="C_A > 0")
        ax.plot(ts[neg], np.abs(CA[neg]), color="tab:red", marker="o", ms=8, lw=0, mfc="none", label="C_A < 0")
        ax.set_yscale("log")
        ax.legend(fontsize=9)
        scaleword = "log"
    else:
        ax.axhline(0.0, color="gray", lw=1, alpha=0.7)
        ax.plot(ts[1:], CA[1:], color="tab:red", marker="o", ms=5, lw=1.2)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        scaleword = "linear"
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$|C_A(t)|$" if LOG else r"$C_A(t)$")
    ax.set_title(r"FREE %s diagram-A $\sigma^2_{00}$ correlator (%s)" % (LTAG, scaleword))
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/diagramA_corr_%s_%s_claude.png" % (scaleword, LTAG)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
