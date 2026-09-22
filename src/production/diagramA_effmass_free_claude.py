#!/usr/bin/env python3
# diagramA_effmass_free_claude.py
#   Effective mass of DIAGRAM A alone (the same-timeslice-contraction Wick term of the sigma^2_00
#   two-meson correlator; Fig 1 of qed3int_v3-4.pdf p.21).  Diagram A = the [4]-cycle whose two sink
#   vertices are cyclically adjacent -> the same-slice equal-time loop = tilde_tau, which turns the pair
#   into a NONLOCAL single-meson interpolator.  Free limit (ENS=free, single exact config): its effmass
#   should plateau at the (2,2) single meson 2E_1 (Delta=4), NOT at the two-meson 2 m_sigma=0.756.
#
#   Reuses the validated per-diagram builder fs_gevp_point_perclass_claude (folded kernels, S+Stilde legs).
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 diagramA_effmass_free_claude.py

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
SPLIT = int(os.environ.get("SPLIT", "1"))


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    nsite = dual.shape[0]

    # single channel: sigma^2_00 (the Kw = Y00-area outer-product kernel), coincident (OFF=(0,0))
    Kw = np.outer(dual * dc.Y00, dual * dc.Y00)
    KER = [Kw]
    OFF = [(0, 0)]

    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s FREE L=%d  DIAGRAM-A effmass of sigma^2_00  (m_sig=%.3f 2m_sig=%.3f  (2,2)~2E1~0.53)"
          % (tag, dc.L, msig, 2 * msig))

    # per-class diagonal correlator for the single config k=0 ; C[class, op, op, dt]
    pc.DTMAX = DTMAX
    C = pc.perclass_one_config(dc.KS[0], KER, OFF, dual)
    iA = pc.CLASSES.index("A")
    cA = C[iA, 0, 0, :].real                                  # diagram-A sigma^2_00 diagonal correlator (contact IN)

    # contact-subtracted diagram A = <O_A O_A>, O_A = psibar tilde_tau psi (tilde_tau = tau_eq - 1/2 I),
    # which is orthogonal to the Delta=2 ground and exposes the (2,2) = {1,1,1,1} at 2E_1.
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    cOA = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        ns = twin - dt
        if ns <= 0:
            continue
        acc = 0.0
        for s in range(ns):
            t = s + dt
            acc += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
        cOA[dt] = acc / ns

    # effective mass, log ratio  m_eff(dt) = log( C(dt)/C(dt+1) )
    dts = np.arange(DTMAX)
    with np.errstate(all="ignore"):
        meff = np.log(cA[:-1] / cA[1:])
        meffOA = np.log(cOA[:-1] / cOA[1:])

    print("\n#  dt |   C_A(contact IN)   m_eff_A |   C(O_A, contact OUT)  m_eff_OA")
    for dt in range(1, DTMAX - 1):
        mm = meff[dt] if np.isfinite(meff[dt]) else np.nan
        mo = meffOA[dt] if np.isfinite(meffOA[dt]) else np.nan
        print("#  %2d | %+.6e   %7.4f |   %+.6e   %7.4f" % (dt, cA[dt], mm, cOA[dt], mo))

    # also print the OTHER classes' effmass at a mid dt for context
    print("\n# other diagrams' m_eff at dt=%d (for contrast):" % (DTMAX // 2))
    for ic, cl in enumerate(pc.CLASSES):
        cc = C[ic, 0, 0, :].real
        with np.errstate(all="ignore"):
            m = np.log(cc[:-1] / cc[1:])
        d = DTMAX // 2
        print("#   diagram %s : m_eff=%7.4f   C=%+.3e" % (cl, m[d] if np.isfinite(m[d]) else np.nan, cc[d]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    # NOTE: both curves have the equal-time contact removed (AblkS subtracts -1/2 I).  The difference is
    # operator structure: diagram A of the 4-fermion sigma^2 retains Delta=2 ground overlap -> m_sigma;
    # the 2-fermion single bilinear O_A (sigma_3-odd) is orthogonal to the ground -> (2,2)=2E_1.
    g = np.isfinite(meff) & (cA[:-1] * cA[1:] > 0)
    ax.plot(dts[:-1][g], meff[g], color="tab:red", marker="o", ms=5, lw=1.2,
            label=r"diagram A of $\sigma^2$ (4-fermion) $\to m_\sigma$ ground")
    go = np.isfinite(meffOA) & (cOA[:-1] * cOA[1:] > 0)
    ax.plot(dts[:-1][go], meffOA[go], color="tab:purple", marker="D", ms=5, lw=1.2,
            label=r"$O_A=\bar\psi\tilde\tau\psi$ (single bilinear) $\to (2,2)$")
    ax.axhline(2 * msig, color="tab:blue", ls="--", lw=1, alpha=0.7)
    ax.text(DTMAX * 0.55, 2 * msig + 0.01, r"$2m_\sigma=%.3f$ (two-meson)" % (2 * msig), fontsize=9, color="tab:blue")
    ax.axhline(0.53, color="k", ls=":", lw=1, alpha=0.6)
    ax.text(DTMAX * 0.55, 0.53 + 0.01, r"$2E_1\approx0.53$  $(2,2)$ single meson", fontsize=9, color="k")
    ax.axhline(msig, color="tab:green", ls=":", lw=1, alpha=0.5)
    ax.text(DTMAX * 0.55, msig + 0.01, r"$m_\sigma=%.3f$" % msig, fontsize=9, color="tab:green")
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(0, DTMAX - 2)
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  DIAGRAM A effmass of $\sigma^2_{00}$ (single-meson channel)" % dc.L)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/diagramA_effmass_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
