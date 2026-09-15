#!/usr/bin/env python3
# gevp_OA_sig22_free_claude.py  [FREE -- {1, O_A, sigma_00^2 (PS^2), O_22} GEVP]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 gevp_OA_sig22_free_claude.py  (or LREF=2 NVDIR=distill_Nv84)
#
# Basis: identity (vacuum) + O_A (psibar tilde_tau psi, sigma3-MIXED single-meson (2,2), mixes with sigma^2)
#        + PS^2 = sigma_00^2 (two-meson) + O_22 = psibar Q2 tilde_tau psi (SHELL-projected (2,2), COUPLING).
# Q2 = V^dag Pi_{lambda2} V is the physical lambda=2 spectral projector from the perambulator time-
# dependence (sigma22_lattice_claude.py): K(dt)=<tau(s+dt,s)>_s = sum_shell Q_shell exp(-E dt); the 2nd
# energy cluster (degeneracy 8 = j=3/2) is lambda=2.  The BARE shell projector psibar Q2 psi is sigma3-EVEN
# and DECOUPLES from sigma^2 (Tr[Gamma GGG]=0, useless); a LOCAL spin (sigma3) insertion does not fix it
# (sigma22_spininsert_test_claude.py).  Combining with the NON-LOCAL tilde_tau (sigma3 tt sigma3 = -tt-1)
# makes O_22 = psibar Q2 tilde_tau psi sigma3-MIXED so it MIXES with sigma^2 (|C2Q/C2A|~0.85) while still
# sharply projecting lambda=2 -> a clean coupling (2,2).  FULL correlators; identity carries the vacuum.
# Ref: qed3int_v3-4.pdf Eq(5.5).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = int(os.environ.get("T0", "3"))
DT_DECOMP = int(os.environ.get("DT_DECOMP", "3"))


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


def shell_projector(tau, twin, Nv, dt_dec):
    ns = twin - dt_dec
    K = np.zeros((Nv, Nv), complex)
    for s in range(ns):
        K += tau[s + dt_dec, s]
    K /= ns
    mu, R = np.linalg.eig(K)
    Rinv = np.linalg.inv(R)
    E = -np.log(np.abs(mu)) / dt_dec
    order = np.argsort(E)
    E = E[order]
    R = R[:, order]
    Rinv = Rinv[order, :]
    # cluster by E; lambda=2 = second cluster
    clusters = []
    i = 0
    while i < Nv:
        j = i
        while j + 1 < Nv and abs(E[j + 1] - E[i]) < 0.02:
            j += 1
        clusters.append((i, j))
        i = j + 1
    i2, j2 = clusters[1]
    sel = list(range(i2, j2 + 1))
    Q2 = R[:, sel] @ Rinv[sel, :]
    return Q2, clusters, E


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  {1, O_A, sigma^2, O_22} GEVP  m_sig~%.3f 2m_sig~%.3f (2,2)~0.52  T0=%d"
          % (tag, dc.L, msig, 2 * msig, T0))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]

    Q2, clusters, Esh = shell_projector(tau, twin, Nv, DT_DECOMP)
    print("# shell degeneracies (should be 4,8,12 = lambda=1,2,3): %s"
          % [c[1] - c[0] + 1 for c in clusters])
    PQ = [Q2 @ tt[a, a] for a in range(twin)]                        # O_22 vertex = Q2 @ tilde_tau (coupling)

    # one-points
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    oQ = np.mean([(-np.trace(PQ[a] @ tt[a, a])).real for a in range(twin)])
    print("# one-points: <sigma^2>=%.4e  <O_A>=%.4e  <O_22>=%.4e" % (o2, oA, oQ))

    CAA = np.zeros(twin); C22 = np.zeros(twin); CQQ = np.zeros(twin)
    C2A = np.zeros(twin); C2Q = np.zeros(twin); CAQ = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        aAA = a22 = aQQ = a2A = a2Q = aAQ = 0.0
        for s in range(ns):
            t = s + dt
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
            aQQ += (-np.trace(PQ[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            aAQ += (-np.trace(PA[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum().real
            a2A += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            a2Q += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PQ[t] @ tau[t, s])).real
        CAA[dt], C22[dt], CQQ[dt] = aAA / ns, a22 / ns, aQQ / ns
        C2A[dt], C2Q[dt], CAQ[dt] = a2A / ns, a2Q / ns, aAQ / ns

    # coupling check: <O_22 sigma^2> connected should now be nonzero (~0.85 |C2A|) since O_22 is sigma3-mixed
    print("# coupling check: |C2Q|/|C2A| at t0 = %.2e / %.2e = %.3f (should be O(1), coupling ON)"
          % (abs(C2Q[T0]), abs(C2A[T0]), abs(C2Q[T0]) / (abs(C2A[T0]) + 1e-300)))

    # FULL 4x4 correlator; order {0:1, 1:O_A, 2:sigma^2, 3:O_22}
    Cts = np.zeros((twin, 4, 4))
    for dt in range(twin):
        Cts[dt] = np.array([
            [1.0,   oA,               o2,               oQ],
            [oA,    CAA[dt] + oA ** 2, C2A[dt] + o2 * oA, CAQ[dt] + oA * oQ],
            [o2,    C2A[dt] + o2 * oA, C22[dt],           C2Q[dt] + o2 * oQ],
            [oQ,    CAQ[dt] + oA * oQ, C2Q[dt] + o2 * oQ, CQQ[dt] + oQ ** 2]])

    lam, nlev = solve_gevp(Cts, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# levels kept = %d" % nlev)
    print("\n  t  |  m0(vac)   m1        m2        m3     [0; (2,2)~0.52; 2m_sig=%.3f]" % (2 * msig))
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
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
        ax.text(twin - 3, y + 0.01, lab, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  $\{1,\ O_A,\ \sigma_{00}^2,\ O_{22}{=}\bar\psi Q_2\tilde\tau\psi\}$ GEVP levels" % dc.L)
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/gevp_OA_sig22tt_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
