#!/usr/bin/env python3
# two_meson_ps2_to_fs_free_claude.py  [FREE -- PS-PS (two-meson) -> FS (single) transition amplitude]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_ps2_to_fs_free_claude.py
#
# The pure-PS 2->1 triangle <sigma_PS(t) sigma_PS^2(0)> ~ 0 (checked).  Here we ask whether the FURNISHED
# scalar sigma_FS = eta^dag xi - xi^dag (1-D_ov^dag) eta overlaps the two-PS-meson state:
#   amplitude  <sigma_FS(t) sigma_PS^2(0)>  (source sigma_PS^2 legs = tau; sink FS vertex furnished by Gamma).
# Kernels as used now: PS leg = tau (= D^{-1}); FS furnished leg = -taugw (= -(1-D_ov^dag)D^{-1}, coincides
# with (1-D_ov^dag)D^{-dag} at L1).  Codebase defs: jj_local_ylm_scalar_conn_stoch_claude.cu (V++, V--^FS).
# Triangle (contraction_v3.pdf) with the sink single-meson return leg swapped tau -> -taugw (two variants
# for which side carries the furnishing -- the thing to confirm).  tilde_tau contact on the source sigma^2.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  PS-PS -> FS transition (kernels: PS=tau, FS=-taugw)  m_sigma~%.3f  2m_sigma~%.3f"
          % (tag, dc.L, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    gw = -taugw                                        # FS leg (channel convention: leg = -taugw)

    triPP = np.zeros(twin)                             # PS^2 -> PS  (baseline, ~0)
    triFSa = np.zeros(twin)                            # PS^2 -> FS, furnish RETURN leg (t,s)
    triFSb = np.zeros(twin)                            # PS^2 -> FS, furnish INCOMING leg (s,t)
    csFSa = np.zeros(twin)                             # FS single two-point ref: -Tr[Phi tau Phi (-taugw)]
    for dt in range(twin):
        ns = twin - dt
        aPP = aFa = aFb = acs = 0.0
        for s in range(ns):
            t = s + dt
            aPP += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s])).real
            aFa += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ gw[t, s])).real
            aFb += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ gw[s, t] @ Phi[t] @ tau[t, s])).real
            acs += (-np.trace(Phi[s] @ tau[s, t] @ Phi[t] @ gw[t, s])).real   # PS-FS single two-point
        triPP[dt], triFSa[dt], triFSb[dt], csFSa[dt] = aPP / ns, aFa / ns, aFb / ns, acs / ns

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))

    ePP, eFa, eFb, ecs = eff(triPP), eff(triFSa), eff(triFSb), eff(csFSa)
    print("\n  dt |  triPP(PS)   triFS_a      triFS_b    |  PS-FS single | m: FSa   FSb   PS-FSsingle")
    for dt in range(1, min(24, twin - 1)):
        print("  %2d | %+.3e  %+.3e  %+.3e | %+.3e |  %5.2f  %5.2f  %5.2f"
              % (dt, triPP[dt], triFSa[dt], triFSb[dt], csFSa[dt], eFa[dt], eFb[dt], ecs[dt]))
    print("# |triPP| max=%.2e (baseline ~0) ; |triFS_a| max=%.2e ; |triFS_b| max=%.2e ; |PS-FS single| max=%.2e"
          % (np.abs(triPP).max(), np.abs(triFSa).max(), np.abs(triFSb).max(), np.abs(csFSa).max()))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(triPP, "tab:gray", "o", r"PS$^2\to$PS (baseline)"),
                            (triFSa, "tab:red", "s", r"PS$^2\to$FS (furnish t,s)"),
                            (triFSb, "tab:blue", "^", r"PS$^2\to$FS (furnish s,t)"),
                            (csFSa, "tab:green", "D", r"PS-FS single 2pt")]:
        ax.semilogy(dts, np.abs(C[1:]), color=col, marker=mk, ms=3, lw=1, label=lab)
    ax.set_xlabel("dt"); ax.set_ylabel(r"$|C(dt)|$")
    ax.set_title(r"FREE L=%d  PS$^2\to$FS transition (kernels tau, -taugw)" % dc.L)
    ax.legend(fontsize=8)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_ps2_to_fs_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
