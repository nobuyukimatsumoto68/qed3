#!/usr/bin/env python3
# two_meson_gevp_OA_free_claude.py  [FREE -- {1, sigma_00, sigma_00^2, O_A} GEVP]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_gevp_OA_free_claude.py
#
# O_A = psibar tilde_tau psi  (non-local scalar; mode vertex Phi_A = Phi_00 @ tilde_tau) overlaps the
# (2,2) excited meson (~0.52) that contaminates sigma_00^2 via diagram A, is orthogonal to sigma_00
# (<O_A sigma_00>=0), and MIXES with sigma_00^2 (<O_A sigma^2> != 0, unlike sigma_00 whose C12=0).
# So the GEVP {1, sigma_00, sigma_00^2, O_A} can project (2,2) out and expose the 2m_sigma two-meson.
# Identity carries the vacuum; FULL correlators (disconnected pieces restored).  Ref: qed3int_v3-4.pdf
# Eq (5.5) two terms; qed3_v2-6.pdf App C free spectrum.

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
    print("# ENS=%s  FREE L=%d  {1, sigma_00, sigma_00^2, O_A} GEVP  m_sig~%.3f 2m_sig~%.3f  (2,2)~0.52"
          % (tag, dc.L, msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]                    # O_A vertex = Phi_00 @ tilde_tau

    # one-points
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])  # <sig^2> (D'_S)
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])                        # <O_A>
    print("# one-points: <sigma_00^2> = %.4e   <O_A> = %.4e" % (o2, oA))

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

    # FULL correlators: restore disconnected <X><Y>.  <sigma_00>=0 so C11,C12,C1A already full.
    # C22 (diags) already full.  Add o2*oA to C2A and oA^2 to CAA.
    Cts = np.zeros((twin, 4, 4))
    for dt in range(twin):
        Cts[dt] = np.array([[1.0,  0.0,      o2,          oA],
                            [0.0,  C11[dt],  C12[dt],     C1A[dt]],
                            [o2,   C12[dt],  C22[dt],     C2A[dt] + o2 * oA],
                            [oA,   C1A[dt],  C2A[dt] + o2 * oA, CAA[dt] + oA ** 2]])

    lam, nlev = solve_gevp(Cts, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# levels kept = %d" % nlev)
    print("\n  t  | m0(vac)   m1(sigma)   m2        m3     [0; 0.378; 0.52(2,2); 0.756(2mes)]")
    for t in range(T0 + 1, twin - 2):
        r = "  ".join("%7.4f" % em[t, i] if i < nlev and np.isfinite(em[t, i]) else "   --- " for i in range(4))
        print("  %2d | %s" % (t, r))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:green", "tab:blue"]
    labs = ["level 0 (vac)", "level 1", "level 2", "level 3"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % 4], marker="o", ms=3, lw=1, label=labs[i])
    for y, lab in [(0.0, "vac"), (msig, r"$m_\sigma$"), (0.52, r"$(2,2)$"), (2 * msig, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.5)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d {$1,\sigma_{00},\sigma_{00}^2,O_A$} GEVP" % dc.L)
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_OA_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
