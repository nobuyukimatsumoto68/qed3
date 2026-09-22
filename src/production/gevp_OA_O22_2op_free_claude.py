#!/usr/bin/env python3
# gevp_OA_O22_2op_free_claude.py  [FREE -- O_A alone, O_22 alone, and {O_A,O_22} 2x2 GEVP]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 gevp_OA_O22_2op_free_claude.py
#
# Purpose (diagnostic for the spurious-level reading):
#   (i)   O_A = psibar tilde_tau psi        single-operator effective mass  (CAA)
#   (ii)  O_22 = psibar Q2 tilde_tau psi    single-operator effective mass  (CQQ)
#   (iii) {O_A, O_22} 2x2 GEVP on the CONNECTED (vacuum-subtracted) matrix
#            C = [[CAA, CAQ], [CAQ, CQQ]]
#         -> level-0 effmass + eigenvector (expect genuine (2,2)=2E1~0.52, eigenvector ~ all O_22),
#            level-1 effmass (expect the spurious direction that decays to null).
# Vertices/contractions identical to gevp_OA_sig22_free_claude.py; only the basis is {O_A,O_22}.
# tilde_tau = tau with contact-subtracted diagonal (GW contact 1/2).  Q2 = V^dag Pi_{lambda2} V.

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


def gevp2(Cts, t0, tol=1e-11):
    # returns eigenvalues sorted desc, eigenvectors in the operator basis, and nlev
    twin = Cts.shape[0]
    C0 = 0.5 * (Cts[t0] + Cts[t0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > tol * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    nlev = int(keep.sum())
    lam = np.full((twin, nlev), np.nan)
    vecs = [None] * twin
    for t in range(twin):
        M = Uk.T @ (0.5 * (Cts[t] + Cts[t].T)) @ Uk
        val, vec = np.linalg.eig(M)
        order = np.argsort(val.real)[::-1]
        lam[t] = val.real[order]
        vecs[t] = Uk @ vec[:, order]
    return lam, vecs, nlev


def effmass(C):
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  FREE L=%d  O_A / O_22 single-op + {O_A,O_22} 2x2 GEVP  T0=%d" % (tag, dc.L, T0))
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
    PQ = [Q2 @ tt[a, a] for a in range(twin)]

    CAA = np.zeros(twin)
    CQQ = np.zeros(twin)
    CAQ = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        aAA = 0.0
        aQQ = 0.0
        aAQ = 0.0
        for s in range(ns):
            t = s + dt
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
            aQQ += (-np.trace(PQ[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            aAQ += (-np.trace(PA[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
        CAA[dt] = aAA / ns
        CQQ[dt] = aQQ / ns
        CAQ[dt] = aAQ / ns

    # single-operator effective masses
    emA = effmass(CAA)
    emQ = effmass(CQQ)
    print("\n# single-operator correlators and effective masses")
    print("#  t  |   CAA         CQQ         CAQ     |  m_A(O_A)   m_Q(O_22)")
    for t in range(1, min(twin - 1, 24)):
        sA = "%8.4f" % emA[t] if np.isfinite(emA[t]) and CAA[t] * CAA[t + 1] > 0 else "  sign/-- "
        sQ = "%8.4f" % emQ[t] if np.isfinite(emQ[t]) and CQQ[t] * CQQ[t + 1] > 0 else "  sign/-- "
        print("#  %2d | %10.3e  %10.3e  %10.3e |  %s  %s" % (t, CAA[t], CQQ[t], CAQ[t], sA, sQ))

    # {O_A, O_22} 2x2 GEVP on the connected (vacuum-subtracted) matrix; order {0:O_A, 1:O_22}
    Cts = np.zeros((twin, 2, 2))
    for dt in range(twin):
        Cts[dt] = np.array([[CAA[dt], CAQ[dt]], [CAQ[dt], CQQ[dt]]])
    lam, vecs, nlev = gevp2(Cts, T0)
    emg = effmass(lam)
    print("\n# {O_A,O_22} 2x2 GEVP (connected)  levels kept = %d" % nlev)
    print("#  t  |  m_lvl0     m_lvl1   |  lvl0 eigvec [O_A, O_22] (normalized)")
    for t in range(T0 + 1, twin - 2):
        r = []
        for i in range(nlev):
            ok = np.isfinite(emg[t, i]) and (lam[t, i] * lam[t + 1, i] > 0)
            r.append("%8.4f" % emg[t, i] if ok else "  ---   ")
        v0 = vecs[t][:, 0]
        v0 = v0 / (np.linalg.norm(v0) + 1e-300)
        print("#  %2d | %s | [%+.3f, %+.3f]" % (t, "  ".join(r), v0[0].real, v0[1].real))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    gA = np.isfinite(emA) & (CAA[:-1] * CAA[1:] > 0)
    gQ = np.isfinite(emQ) & (CQQ[:-1] * CQQ[1:] > 0)
    ax.plot(ts[gA], emA[gA], color="tab:red", marker="o", ms=4, lw=1, label=r"$O_A$ alone")
    ax.plot(ts[gQ], emQ[gQ], color="tab:blue", marker="s", ms=4, lw=1, label=r"$O_{22}$ alone")
    for i in range(nlev):
        g = np.isfinite(emg[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        col = ["k", "tab:gray"][i % 2]
        mk = ["^", "v"][i % 2]
        ax.plot(ts[g], emg[g, i], color=col, marker=mk, ms=3, lw=1, ls="--",
                label=r"GEVP level %d" % i)
    for y, lab in [(0.378, r"$2m_\sigma^{(1shell)}$"), (0.52, r"$(2,2){=}2E_1$"), (0.756, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls=":", lw=0.8, alpha=0.4)
        ax.text(twin - 3, y + 0.008, lab, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  $O_A$ / $O_{22}$ single-op + $\{O_A,O_{22}\}$ GEVP" % dc.L)
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/gevp_OA_O22_2op_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
