#!/usr/bin/env python3
# point_split_test_claude.py  -- test the covariant point-split scalar O_ps.
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 point_split_test_claude.py
# Checks: (a) mixing <O_ps(t) sigma_00^2(0)>  (must be != 0, unlike sigma_00 whose C12=0),
#         (b) 2-point <O_ps O_ps> effmass (which single-meson states; is (2,2)~0.52 there),
#         (c) <O_ps sigma_00> single mixing.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import geom_hopping_claude as gh


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  point-split O_ps test  m_sig~%.3f 2m_sig~%.3f  (2,2)~0.52"
          % (tag, dc.L, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = tau.shape[-1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    # vertices
    Wps, nns = gh.pointsplit_vertex(dc.GEOM, dc.L, dual, dc.Y00, dc.NS)
    Phi = []
    Pps = []
    for a in range(twin):
        Vt = V[tsrc0 + a].T                              # (2Ns, Nv)
        Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        Pps.append(Vt.conj().T @ Wps @ Vt)
    # <O_ps> one-point (should ~0: no scalar condensate at a link, parity)
    o_ps = np.mean([(-np.trace(Pps[a] @ tt[a, a])).real for a in range(twin)])
    print("# <O_ps> one-point = %.4e" % o_ps)

    mix = np.zeros(twin)     # <O_ps sigma^2>
    Cpp = np.zeros(twin)     # <O_ps O_ps>
    Cps0 = np.zeros(twin)    # <O_ps sigma_00>
    for dt in range(twin):
        ns = twin - dt
        am = ac = a0 = 0.0
        for s in range(ns):
            t = s + dt
            am += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Pps[t] @ tau[t, s])).real
            ac += (-np.trace(Pps[t] @ tau[t, s] @ Pps[s] @ tau[s, t])).real
            a0 += (-np.trace(Pps[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
        mix[dt], Cpp[dt], Cps0[dt] = am / ns, ac / ns, a0 / ns

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))

    emix, epp, e0 = eff(mix), eff(Cpp), eff(Cps0)
    print("\n  dt |  <O_ps sig^2>(mix)  m_mix | <O_ps O_ps>   m_pp | <O_ps sig00>  m_0")
    for dt in range(1, min(22, twin - 1)):
        print("  %2d |  %+.3e      %5.2f | %+.3e   %5.2f | %+.3e  %5.2f"
              % (dt, mix[dt], emix[dt], Cpp[dt], epp[dt], Cps0[dt], e0[dt]))
    print("# |<O_ps sigma^2>| max = %.3e   (nonzero => usable to project a state via GEVP)" % np.abs(mix).max())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(Cpp, "tab:red", "o", r"$\langle O_\mathrm{ps}O_\mathrm{ps}\rangle$"),
                            (mix, "tab:blue", "s", r"$\langle O_\mathrm{ps}\sigma^2\rangle$ mixing"),
                            (Cps0, "tab:green", "^", r"$\langle O_\mathrm{ps}\sigma_{00}\rangle$")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=3, lw=1, label=lab)
    for y, lab in [(msig, r"$m_\sigma=0.38$"), (0.52, r"$(2,2)=0.52$"), (2 * msig, r"$2m_\sigma=0.76$")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("dt"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d point-split $O_\mathrm{ps}$: 2pt, mixing" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/point_split_test_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
