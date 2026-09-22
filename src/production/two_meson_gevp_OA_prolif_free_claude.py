#!/usr/bin/env python3
# two_meson_gevp_OA_prolif_free_claude.py  [FREE -- {1, sigma_00^2, O_A} proliferated, step DT_SHIFT]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 DT_SHIFT=2 python3 two_meson_gevp_OA_prolif_free_claude.py
#
# sigma_00 dropped (it is fully decoupled: C12=C1A=<sigma>=0).  Two-meson basis = {1, sigma_00^2, O_A}
# with O_A = psibar tilde_tau psi (mode vertex Phi_00 @ tilde_tau).  Proliferate the non-identity ops with
# a time shift of DT_SHIFT steps (asymm_gevp.pdf; Wagman 2406.20009): {sigma^2,O_A, sigma^2^(D),O_A^(D)} 5x5,
# block-Hankel c(t), c(t+D), c(t+2D).  A larger shift D gives more-independent inflated ops -> steadier plateau.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = int(os.environ.get("T0", "3"))
D = int(os.environ.get("DT_SHIFT", "2"))
NLEV = int(os.environ.get("NLEV", "0"))         # >0: rebase (truncate C(t0) to its top NLEV eigenvalues)
NOVAC = int(os.environ.get("NOVAC", "0"))       # 1: drop the identity; use vacuum-subtracted (connected) c


def solve_gevp(Cts, T0, tmax, tol=1e-11):
    n = Cts.shape[1]
    C0 = 0.5 * (Cts[T0] + Cts[T0].T)
    wv, Uv = np.linalg.eigh(C0)
    order = np.argsort(wv)[::-1]                 # descending eigenvalue
    if NLEV > 0:
        sel = order[:NLEV]                        # top NLEV (rebase into NLEV states)
    else:
        sel = order[wv[order] > tol * wv.max()]
    Uk = Uv[:, sel] / np.sqrt(wv[sel])
    nlev = len(sel)
    lam = np.full((Cts.shape[0], nlev), np.nan)
    for t in range(tmax):
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
    print("# ENS=%s  FREE L=%d  {1, sigma_00^2, O_A} proliferated  shift D=%d  m_sig~%.3f 2m_sig~%.3f (2,2)~0.52"
          % (tag, dc.L, D, msig, 2 * msig))
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
    o = np.array([o2, oA])

    C22 = np.zeros(twin); C2A = np.zeros(twin); CAA = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a22 = a2A = aAA = 0.0
        for s in range(ns):
            t = s + dt
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum().real
            a2A += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
        C22[dt], C2A[dt], CAA[dt] = a22 / ns, a2A / ns, aAA / ns

    # base 2x2 c(t) for {sigma^2, O_A} (FULL; disconnected restored)
    c = np.zeros((twin, 2, 2))
    for t in range(twin):
        c[t] = np.array([[C22[t],            C2A[t] + o2 * oA],
                         [C2A[t] + o2 * oA,  CAA[t] + oA ** 2]])
    if NOVAC:
        oo = np.outer(o, o)
        c = c - oo[None, :, :]                    # connected (vacuum removed); no identity operator

    # inflated:  NOVAC -> 4x4 block-Hankel of connected c {sigma^2,O_A, sigma^2^(D),O_A^(D)};
    #            else   -> 5x5 with the identity carrying the vacuum.
    def M(t):
        if NOVAC:
            m = np.zeros((4, 4))
            m[0:2, 0:2] = c[t]
            m[0:2, 2:4] = c[t + D]
            m[2:4, 0:2] = c[t + D]
            m[2:4, 2:4] = c[t + 2 * D]
            return m
        m = np.zeros((5, 5))
        m[0, 0] = 1.0
        for i in range(2):
            m[0, 1 + i] = m[1 + i, 0] = o[i]
            m[0, 3 + i] = m[3 + i, 0] = o[i]
        m[1:3, 1:3] = c[t]
        m[1:3, 3:5] = c[t + D]
        m[3:5, 1:3] = c[t + D]
        m[3:5, 3:5] = c[t + 2 * D]
        return m

    tmax = twin - 2 * D
    ndim = 4 if NOVAC else 5
    Cts = np.array([M(t) if t < tmax else np.zeros((ndim, ndim)) for t in range(twin)])
    lam, nlev = solve_gevp(Cts, T0, tmax)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# levels kept = %d" % nlev)
    print("\n  t  |  " + "  ".join("m%d" % i for i in range(nlev)) + "    [0; 0.378; 0.52; 0.756]")
    for t in range(T0 + 1, tmax - 1):
        row = "  ".join("%6.3f" % em[t, i] if i < nlev and np.isfinite(em[t, i]) else "  --- " for i in range(nlev))
        print("  %2d |  %s" % (t, row))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, tmax - 1)
    cols = ["tab:gray", "tab:green", "tab:blue", "tab:purple", "tab:orange"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:tmax - 1, i]) & (lam[:tmax - 1, i] * lam[1:tmax, i] > 0)
        ax.plot(ts[g], em[:tmax - 1][g, i], color=cols[i % len(cols)], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, msig, 0.52, 2 * msig]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d {$1,\sigma^2,O_A$} proliferated (shift D=%d)" % (dc.L, D))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_OA_prolif_D%d_t0%d_free_L%d_claude.png" % (D, T0, dc.L)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
