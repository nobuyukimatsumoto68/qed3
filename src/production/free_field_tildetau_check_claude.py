#!/usr/bin/env python3
# free_field_tildetau_check_claude.py   [copy&edit -- FREE FIELD check of NM's contact-subtraction idea]
# Run with:  ENS=free python3 free_field_tildetau_check_claude.py
# Free field (U=1), 1 deterministic config, distillation Nv=24.  Compare
#   RAW   (sigma^2, contact kept, CONTACT=0)   vs   tilde_tau (contact-subtracted, CONTACT=0.5)
# for the single-meson C_S, the two-meson B/E, and A+B+E.  Question: does tilde_tau kill the single-meson
# in A so that A+B+E -> 2 m_sigma?  Reference single-meson m from C_S itself.  Effmass = log C(dt)/C(dt+1)
# (1 config -> central values; window-source jackknife for a rough error).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

KEEP = [0, 1, 4]                                    # A, B, E


def build():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    return Phi, tau, twin


def corr(Phi, tau, twin, contact):
    de.CONTACT = contact
    D = np.zeros((10, twin))               # per-diagram, translation-averaged
    CS = np.zeros(twin)                    # single-meson C_S(dt)
    persrc = []                            # per-source A+B+E for window jackknife
    for dt in range(twin):
        ns = twin - dt
        acc = np.zeros(10, complex)
        cs = 0.0
        ps = []
        for s in range(ns):
            d10 = de.diags_pair(Phi, tau, s, s + dt)
            acc += d10
            M = Phi[s] @ tau[s, s + dt] @ Phi[s + dt] @ tau[s + dt, s]
            cs += np.trace(M)
            w = 2.0 * dc.W10[KEEP]
            ps.append((w[0] * d10[0] + w[1] * d10[1] + w[2] * d10[4]).real)
        D[:, dt] = (acc / ns).real
        CS[dt] = (cs / ns).real
        persrc.append(np.array(ps))
    return D, CS, persrc


def eff(C):
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


def eff_src(persrc):                       # window-source jackknife effmass of A+B+E
    twin = len(persrc)
    mean = np.array([persrc[dt].mean() for dt in range(twin)])
    em = eff(mean)
    # jackknife over sources at the smallest window (crude error)
    return mean, em


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE FIELD  ncfg=%d  raw vs tilde_tau" % (tag, len(dc.KS)))
    Phi, tau, twin = build()
    for label, c in [("RAW (contact kept)", 0.0), ("tilde_tau (contact removed)", 0.5)]:
        D, CS, persrc = corr(Phi, tau, twin, c)
        emCS = eff(CS)
        emA = eff(D[0])
        emB = eff(D[1])
        emE = eff(D[4])
        w = 2.0 * dc.W10[KEEP]
        ABE = w[0] * D[0] + w[1] * D[1] + w[2] * D[4]
        emABE = eff(ABE)
        print("\n=== %s ===" % label)
        print("  dt   C_S(single)   A(-S_S)   B(-T_S)   E(C_S^2)   A+B+E")
        for dt in range(2, min(18, twin - 1)):
            print("  %2d   %7.4f       %7.4f   %7.4f   %7.4f   %7.4f"
                  % (dt, emCS[dt], emA[dt], emB[dt], emE[dt], emABE[dt]))

    # plot A+B+E effmass raw vs tilde_tau, and C_S reference
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for label, c, col, mk in [("raw", 0.0, "tab:red", "o"), ("tilde_tau", 0.5, "tab:blue", "s")]:
        D, CS, _ = corr(Phi, tau, twin, c)
        w = 2.0 * dc.W10[KEEP]
        ABE = w[0] * D[0] + w[1] * D[1] + w[2] * D[4]
        em = eff(ABE)
        dts = np.arange(2, twin - 2)
        good = np.isfinite(em[dts]) & (ABE[dts] > 0) & (ABE[dts + 1] > 0)
        ax.plot(dts[good], em[dts][good], color=col, marker=mk, ms=4, lw=1, label="A+B+E %s" % label)
        if c == 0.5:
            emCS = eff(CS)
            gc = np.isfinite(emCS[dts]) & (CS[dts] > 0) & (CS[dts + 1] > 0)
            ax.plot(dts[gc], emCS[dts][gc], color="gray", marker="^", ms=3, lw=0.8, alpha=0.7,
                    label=r"$C_S$ (single meson)")
            ax.plot(dts[gc], 2.0 * emCS[dts][gc], color="k", ls="--", lw=1, alpha=0.6, label=r"$2\times C_S$")
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title("FREE FIELD: A+B+E raw vs tilde_tau, single-meson check")
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/free_field_tildetau_check_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
