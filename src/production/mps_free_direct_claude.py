#!/usr/bin/env python3
# mps_free_direct_claude.py
#   Direct free-L1 single-meson (ell=0, Gamma=1) two-point effmass = m_PS.
#   C_S(dt) = < -Tr[ Phi(t) tau(t,s) Phi(s) tau(s,t) ] >_s  ,  Phi = V^dag diag(A Y00) V (ell=0 wall).
#   Expected: plateau at 2E0 = m_PS = 0.378 (free L1).  Contact IRRELEVANT here (off-diagonal legs only).
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 mps_free_direct_claude.py

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


def main():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    C = np.full(twin, np.nan)
    for dt in range(twin):
        ns = twin - dt
        if ns <= 0:
            continue
        acc = 0.0
        for s in range(ns):
            t = s + dt
            acc += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
        C[dt] = acc / ns

    with np.errstate(all="ignore"):
        meff = np.log(C[:-1] / C[1:])

    print("# FREE L=1 single-meson (ell=0, Gamma=1) C_S effmass -- expect m_PS = 2E0 = 0.378")
    print("#  dt |   C_S(dt)        m_eff")
    for dt in range(1, twin - 1):
        if not np.isfinite(C[dt]):
            continue
        mm = meff[dt] if np.isfinite(meff[dt]) else np.nan
        print("#  %2d | %+.6e   %7.4f" % (dt, C[dt], mm))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(twin - 1)
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    ax.axhline(0.378, color="tab:green", ls="--", lw=1, alpha=0.7)
    ax.text(twin * 0.5, 0.383, r"$m_{PS}=2E_0=0.378$", fontsize=10, color="tab:green")
    g = np.isfinite(meff) & (C[:-1] * C[1:] > 0)
    ax.plot(ts[g], meff[g], color="tab:red", marker="o", ms=5, lw=1.2, label=r"$C_S$ single meson")
    ax.set_ylim(0.2, 0.7)
    ax.set_xlim(0, twin - 2)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=1 single-meson $C_S$ effmass $\to m_{PS}$")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/mps_free_direct_claude.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
