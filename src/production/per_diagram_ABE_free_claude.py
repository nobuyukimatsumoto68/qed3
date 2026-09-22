#!/usr/bin/env python3
# per_diagram_ABE_free_claude.py
#   Per-diagram effmass of the sigma^2_00 diagonal for the SURVIVING connected diagrams A, B, E
#   (the tadpole diagrams C/D/G vanish after the -1/2 GW contact subtraction).  Clean contact-subtracted
#   S-leg (AblkS), single free config.  A = sink-adjacent 4-cycle, B = sink-non-adjacent 4-cycle,
#   E = 2+2 (two mesons).
#   Refs (env, set per L): m_PS, (diagram A)=2E1, two-meson=2 m_PS.
#   Run: ENS=free LREF=2 NVDIR=distill_Nv84 M_PS=0.393 STATE_A=0.69 TWO_MESON=0.786 python3 per_diagram_ABE_free_claude.py

import os
os.environ.setdefault("ENS", "free")
os.environ.setdefault("LREF", "1")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G
import fs_gevp_point_perclass_claude as pc

DTMAX = int(os.environ.get("DTMAX", "24"))
LTAG = os.environ.get("LTAG", "L1")
M_PS = float(os.environ.get("M_PS", "0.378"))
STATE_A = float(os.environ.get("STATE_A", "0.556"))
TWO_MESON = float(os.environ.get("TWO_MESON", "0.756"))
CLS = ["A", "B", "E"]


def main():
    AblkS, AblkSt, twin, nsite, U, tau = G.make_config(dc.KS[0])
    dual = dc.dual_areas_from_mesh().astype(float)
    wY = dual * dc.Y00
    Pmap = G.antipodal_map()
    cls_of = [pc.classify(cy) for cy in G.PERMS]

    C = {c: np.full(DTMAX, np.nan) for c in CLS}
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        vt = [dt, dt, 0, 0]
        offs = set((vt[a], vt[b]) for a in range(4) for b in range(4))
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        vspec = G.op_vspec(0, ('i', 'j'), dual, wY) + G.op_vspec(0, ('k', 'l'), dual, wY)
        acc = {c: 0.0 for c in CLS}
        for ip, cyc in enumerate(G.PERMS):
            cl = cls_of[ip]
            if cl in CLS:
                acc[cl] += G.perm_contrib_folded(cyc, vt, bAS, vspec, Pmap).real
        for c in CLS:
            C[c][dt] = acc[c] / len(s0s)

    eff = {}
    for c in CLS:
        with np.errstate(all="ignore"):
            eff[c] = np.log(C[c][:-1] / C[c][1:])

    print("# %s per-diagram sigma^2_00 effmass (S-leg, contact-subtracted); refs m_PS=%.3f (2,2)=%.3f 2mes=%.3f"
          % (LTAG, M_PS, STATE_A, TWO_MESON))
    print("#  dt |   A         B         E")
    for dt in range(1, DTMAX - 1):
        r = "  ".join("%7.4f" % eff[c][dt] if np.isfinite(eff[c][dt]) else "  ---  " for c in CLS)
        print("#  %2d | %s" % (dt, r))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(DTMAX - 1)
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    for y, lab, col in [(M_PS, r"$m_{PS}=%.3f$" % M_PS, "tab:gray"),
                        (STATE_A, r"(diagram A)$=%.3f$" % STATE_A, "tab:orange"),
                        (TWO_MESON, r"two-meson$=%.3f$" % TWO_MESON, "tab:blue")]:
        ax.axhline(y, color=col, ls="--", lw=1, alpha=0.6)
        ax.text(DTMAX * 0.60, y + 0.008, lab, fontsize=9, color=col)
    cols = {"A": "tab:red", "B": "tab:purple", "E": "tab:green"}
    mkr = {"A": "o", "B": "s", "E": "^"}
    for c in CLS:
        g = np.isfinite(eff[c]) & (C[c][:-1] * C[c][1:] > 0)
        ax.plot(ts[g], eff[c][g], color=cols[c], marker=mkr[c], ms=5, lw=1.1, label="diagram %s" % c)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(0, DTMAX - 2)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE %s per-diagram $\sigma^2_{00}$ effmass: A, B, E" % LTAG)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/per_diagram_ABE_free_%s_claude.png" % LTAG
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
