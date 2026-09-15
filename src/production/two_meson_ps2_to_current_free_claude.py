#!/usr/bin/env python3
# two_meson_ps2_to_current_free_claude.py  [FREE -- PS-PS (two-meson) -> Vector/Axial current ell=0]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_ps2_to_current_free_claude.py
#
# Test whether the ell=0 local current (Pauli sigma^a vertex; jj_local_ylm: tp=s3, sp=(s1+s2)/2) overlaps
# the two-scalar-meson channel -- i.e. whether the ~0.56 "state" seen in diagram A is this (lattice-artifact,
# continuum-forbidden 0++ <-> 1) current.  Triangle  <J^a(t) sigma_PS^2(0)> with the current at the sink:
#   -2 Tr[Phi(s) tilde_tau(s,s) Phi(s) tau(s,t) Phi^a(t) tau(t,s)]  (a=1,2,3).
# Phi^a(t) = V(t)^dag (W_00 (x) sigma^a) V(t)  (Pauli spin vertex, ell=0 wall).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

PAULI = {1: np.array([[0, 1], [1, 0]], complex),
         2: np.array([[0, -1j], [1j, 0]], complex),
         3: np.array([[1, 0], [0, -1]], complex)}


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  PS-PS -> current ell=0 (Pauli s^a)  m_sig~%.3f 2m_sig~%.3f  diagA-state~0.56"
          % (tag, dc.L, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w = dual * dc.Y00                                     # per-site weight (Nsites)
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = tau.shape[-1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    ns_sites = len(dual)

    # scalar wall vertex Phi (identity spin) and current vertices Phi^a (Pauli spin), per timeslice
    Phi = []
    Pa = {1: [], 2: [], 3: []}
    for a in range(twin):
        Vr = V[tsrc0 + a].reshape(Nv, ns_sites, dc.NS)   # [mode, site, spin]
        Vt = V[tsrc0 + a].T
        Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        for aa in range(1, 4):
            WVr = np.einsum('x, st, jxt -> jxs', w, PAULI[aa], Vr, optimize=True)
            Pa[aa].append(np.einsum('kxs, jxs -> kj', np.conj(Vr), WVr, optimize=True))

    tri = {1: np.zeros(twin), 2: np.zeros(twin), 3: np.zeros(twin)}
    for dt in range(twin):
        nsd = twin - dt
        acc = {1: 0.0, 2: 0.0, 3: 0.0}
        for s in range(nsd):
            t = s + dt
            base = Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t]     # source sigma^2 block up to sink
            for aa in range(1, 4):
                acc[aa] += (-2.0 * np.trace(base @ Pa[aa][t] @ tau[t, s])).real
        for aa in range(1, 4):
            tri[aa][dt] = acc[aa] / nsd

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))

    print("\n  dt |  s1-current   s2-current   s3-current(tp) |  m_eff: s1     s2     s3")
    e = {aa: eff(tri[aa]) for aa in range(1, 4)}
    for dt in range(1, min(24, twin - 1)):
        print("  %2d | %+.3e  %+.3e  %+.3e |  %5.2f  %5.2f  %5.2f"
              % (dt, tri[1][dt], tri[2][dt], tri[3][dt], e[1][dt], e[2][dt], e[3][dt]))
    for aa in range(1, 4):
        print("# |PS^2 -> s%d| max = %.3e" % (aa, np.abs(tri[aa]).max()))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for aa, col, mk in [(1, "tab:red", "o"), (2, "tab:green", "^"), (3, "tab:blue", "s")]:
        ee = eff(tri[aa])
        g = np.isfinite(ee) & (tri[aa][:-1] * tri[aa][1:] > 0)
        ax.plot(dts[:len(ee)][g], ee[g], color=col, marker=mk, ms=3, lw=1, label=r"$\sigma^%d$ current" % aa)
    for y, lab in [(msig, r"$m_\sigma$"), (2 * msig, r"$2m_\sigma$"), (0.56, "diagA ~0.56")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("dt"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  PS-PS $\to$ current $\ell=0$ (Pauli $\sigma^a$)" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_ps2_to_current_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
