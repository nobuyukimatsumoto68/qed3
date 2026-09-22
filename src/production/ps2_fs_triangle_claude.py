#!/usr/bin/env python3
# ps2_fs_triangle_claude.py
#   <sigma_FS(t) sigma_PS^2(0)> triangle: does the two-meson sigma_PS^2 couple to a SINGLE sigma_FS?
#   (The PS analogue <sigma_00 sigma^2> = 0 by sigma3-herm; the FS furnishing's sigma3-odd part is NOT
#   protected, so this is expected nonzero -- an O(a) parity-impurity probe.)  Recipe: Fin's
#   fs_ps2_triangle_recipe_claude.md.  Boxed expression (only the src2->sink leg is furnished):
#     <sigma_FS sigma_PS^2> = -Tr[ Phi_0 . ttilde_00 . Phi_0 . tau[s,t] . Phi_t . (-taugw[t,s]) ]   (t=s+dt)
#   ttilde_00 = tau(s,s)-1/2 (plain PS source contact); prefactor -1 (single loop, S3=0).
#   Validation (1): taugw->tau must give 0 (reproduces the protected all-plain triangle).
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 ps2_fs_triangle_claude.py
#        ENS=free LREF=2 NVDIR=distill_Nv84 python3 ps2_fs_triangle_claude.py

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

LTAG = os.environ.get("LTAG", "L%d" % dc.L)


def build():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv

    Css = np.full(twin, np.nan)
    Cfs = np.full(twin, np.nan)          # furnished triangle
    Cval = np.full(twin, np.nan)         # validation: taugw->tau (must be ~0)
    for dt in range(twin):
        ns = twin - dt
        if ns <= 0:
            continue
        ass = 0.0
        afs = 0.0
        aval = 0.0
        for s in range(ns):
            t = s + dt
            ass += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            # boxed: -Tr[ Phi_s ttilde_ss Phi_s tau[s,t] Phi_t (-taugw[t,s]) ]  (sink-arriving leg t<-s furnished)
            afs += (-np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ (-taugw[t, s]))).real
            aval += (-np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ (-tau[t, s]))).real
        Css[dt] = ass / ns
        Cfs[dt] = afs / ns
        Cval[dt] = aval / ns
    return Css, Cfs, Cval, twin


def main():
    Css, Cfs, Cval, twin = build()
    R = Cfs / Css
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    with np.errstate(all="ignore"):
        em = np.log(np.abs(Cfs[:-1] / Cfs[1:]))
    print("# %s <sigma_FS(t) sigma_PS^2(0)> triangle  (m_sig=%.3f)" % (LTAG, msig))
    print("# VALIDATION (taugw->tau, must be ~0 rel. to Css):  max|Cval/Css| = %.2e"
          % np.nanmax(np.abs(Cval / Css)))
    print("#  dt |   C_ss          C_FS(triangle)   R=C_FS/Css    effmass(C_FS)    Cval(->0)")
    for dt in range(1, twin - 1):
        if np.isfinite(Css[dt]) and abs(Css[dt]) > 0:
            e = em[dt] if np.isfinite(em[dt]) else np.nan
            print("#  %2d | %+.5e   %+.5e   %+.6f   %7.4f    %+.2e" % (dt, Css[dt], Cfs[dt], R[dt], e, Cval[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(twin - 1)
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    ax.axhline(2 * msig, color="tab:blue", ls="--", lw=1, alpha=0.6)
    ax.text(twin * 0.5, 2 * msig + 0.01, r"$2E_0=%.3f$" % (2 * msig), fontsize=9, color="tab:blue")
    ax.axhline(msig, color="tab:green", ls=":", lw=1, alpha=0.5)
    ax.text(twin * 0.5, msig + 0.01, r"$m_\sigma=%.3f$" % msig, fontsize=9, color="tab:green")
    g = np.isfinite(em) & (Cfs[:-1] * Cfs[1:] > 0)
    ax.plot(ts[g], em[g], color="tab:red", marker="o", ms=5, lw=1.2, label=r"$\langle\sigma_{FS}\sigma_{PS}^2\rangle$ effmass")
    ax.set_ylim(0.0, 1.2)
    ax.set_xlim(0, twin - 2)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE %s: $\langle\sigma_{FS}(t)\,\sigma_{PS}^2(0)\rangle$ (PS$^2\to$FS single-meson)" % LTAG)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/ps2_fs_triangle_%s_claude.png" % LTAG
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
