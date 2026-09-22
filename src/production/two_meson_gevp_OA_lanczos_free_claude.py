#!/usr/bin/env python3
# two_meson_gevp_OA_lanczos_free_claude.py  [FREE -- {1, sigma_00, sigma_00^2, O_A} + Lanczos inflation]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_gevp_OA_lanczos_free_claude.py
#
# Proliferate the WORKING basis {1, sigma_00, sigma_00^2, O_A} by time-shifting the non-identity operators
# (asymm_gevp.pdf; Wagman 2406.20009): {sigma,sigma^2,O_A, sigma^(1),sigma^2^(1),O_A^(1)} -> 7x7.
# Block-Hankel: base 3x3 c(t) for {sigma,sigma^2,O_A}; shift source or sink by one step -> c(t+1); both -> c(t+2).
# Identity NOT inflated (carries the vacuum).  Goal: extend/steady the level-3 two-meson (2 m_sigma) plateau.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = 3


def solve_gevp(Cts, T0, tol=1e-11):
    twin = Cts.shape[0]
    C0 = 0.5 * (Cts[T0] + Cts[T0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > tol * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    nlev = int(keep.sum())
    lam = np.full((twin, nlev), np.nan)
    for t in range(twin):
        M = Uk.T @ (0.5 * (Cts[t] + Cts[t].T)) @ Uk
        try:
            lam[t] = np.sort(np.linalg.eigvals(M).real)[::-1]
        except Exception:
            pass
    return lam, nlev


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  {1,sigma,sigma^2,O_A}+Lanczos (7x7)  m_sig~%.3f 2m_sig~%.3f (2,2)~0.52"
          % (tag, dc.L, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]

    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    ovec = np.array([0.0, o2, oA])                              # <sigma>, <sigma^2>, <O_A>

    C11 = np.zeros(twin); C12 = np.zeros(twin); C22 = np.zeros(twin)
    C1A = np.zeros(twin); C2A = np.zeros(twin); CAA = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a11 = a12 = a22 = a1A = a2A = aAA = 0.0
        for s in range(ns):
            t = s + dt
            a11 += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            a12 += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Phi[t] @ tau[t, s])).real
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum().real
            a1A += (-np.trace(Phi[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
            a2A += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
        C11[dt], C12[dt], C22[dt] = a11 / ns, a12 / ns, a22 / ns
        C1A[dt], C2A[dt], CAA[dt] = a1A / ns, a2A / ns, aAA / ns

    # FULL 3x3 base correlator c(t) for {sigma, sigma^2, O_A} (disconnected restored)
    c = np.zeros((twin, 3, 3))
    for t in range(twin):
        c[t] = np.array([[C11[t], C12[t],            C1A[t]],
                         [C12[t], C22[t],            C2A[t] + o2 * oA],
                         [C1A[t], C2A[t] + o2 * oA,  CAA[t] + oA ** 2]])

    # 7x7 inflated: [id, sigma, sigma^2, O_A, sigma^(1), sigma^2^(1), O_A^(1)]
    def M(t):
        m = np.zeros((7, 7))
        m[0, 0] = 1.0
        for i in range(3):
            m[0, 1 + i] = m[1 + i, 0] = ovec[i]
            m[0, 4 + i] = m[4 + i, 0] = ovec[i]
        m[1:4, 1:4] = c[t]
        m[1:4, 4:7] = c[t + 1]
        m[4:7, 1:4] = c[t + 1]
        m[4:7, 4:7] = c[t + 2]
        return m

    twin2 = twin - 2
    Cts = np.array([M(t) if t < twin2 else np.zeros((7, 7)) for t in range(twin)])
    lam, nlev = solve_gevp(Cts[:twin2], T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# levels kept = %d" % nlev)
    print("\n  t  |  " + "  ".join("m%d" % i for i in range(nlev)) + "    [0;0.378;0.52;0.756]")
    for t in range(T0 + 1, twin2 - 2):
        row = "  ".join("%6.3f" % em[t, i] if i < nlev and np.isfinite(em[t, i]) else "  --- " for i in range(nlev))
        print("  %2d |  %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin2 - 1)
    cols = ["tab:gray", "tab:red", "tab:green", "tab:blue", "tab:purple", "tab:orange", "tab:brown"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % len(cols)], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, msig, 0.52, 2 * msig]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d {$1,\sigma,\sigma^2,O_A$}+Lanczos 7x7" % dc.L)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_OA_lanczos_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
