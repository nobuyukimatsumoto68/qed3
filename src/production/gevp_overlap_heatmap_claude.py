#!/usr/bin/env python3
# gevp_overlap_heatmap_claude.py
#   Operator->state overlap "wavefunction" matrix for the {1, O_A, sigma^2, O_22^tt} GEVP, at t=20.
#   Solve C(t) v = lam C(t0) v at t0=3; eigenvectors v^n (normalized (v^n)^dag C0 v^n = 1); overlap
#   Z_i^n = (C0 v^n)_i.  Heatmap = |Z_i^n| column-normalized (each STATE's composition over operators).
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 gevp_overlap_heatmap_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de
import sigma22_spininsert_test_claude as st

T0 = int(os.environ.get("T0", "3"))
TMEAS = int(os.environ.get("TMEAS", "20"))


def main():
    de.CONTACT = 0.5
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    Q2 = st.shell_projector(tau, twin, Nv, 3)
    PQ = [Q2 @ tt[a, a] for a in range(twin)]

    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    oQ = np.mean([(-np.trace(PQ[a] @ tt[a, a])).real for a in range(twin)])

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

    def Cmat(dt):
        return np.array([
            [1.0,   oA,               o2,               oQ],
            [oA,    CAA[dt] + oA ** 2, C2A[dt] + o2 * oA, CAQ[dt] + oA * oQ],
            [o2,    C2A[dt] + o2 * oA, C22[dt],           C2Q[dt] + o2 * oQ],
            [oQ,    CAQ[dt] + oA * oQ, C2Q[dt] + o2 * oQ, CQQ[dt] + oQ ** 2]])

    C0 = 0.5 * (Cmat(T0) + Cmat(T0).T)
    Ct = 0.5 * (Cmat(TMEAS) + Cmat(TMEAS).T)
    Ct1 = 0.5 * (Cmat(TMEAS + 1) + Cmat(TMEAS + 1).T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > 1e-11 * wv.max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    nlev = int(keep.sum())
    # generalized eigenproblem in whitened basis at TMEAS
    M = Uk.T @ Ct @ Uk
    lam, u = np.linalg.eigh(M)
    # energies from lam(t)/lam(t+1)
    M1 = Uk.T @ Ct1 @ Uk
    lam1 = np.sort(np.linalg.eigvalsh(M1))[::-1]
    order = np.argsort(lam)[::-1]           # descending lam = ascending energy
    lam = lam[order]
    u = u[:, order]
    with np.errstate(all="ignore"):
        En = np.log(lam / lam1[np.argsort(lam1)[::-1]])
    # overlaps
    ops = ["1", "O_A", r"$\sigma^2$", r"$O_{22}^{\tilde\tau}$"]
    Z = np.zeros((4, nlev))
    for n in range(nlev):
        vn = Uk @ u[:, n]                    # operator-basis eigenvector, (v^n)^dag C0 v^n = 1
        Zn = C0 @ vn                          # overlap Z_i^n
        Z[:, n] = Zn.real
    Zc = Z / np.sqrt((Z ** 2).sum(axis=0, keepdims=True))   # column-normalized composition
    np.set_printoptions(precision=3, suppress=True)
    print("# t0=%d  t=%d   state energies (effmass at t): %s" % (T0, TMEAS, np.array2string(En, precision=3)))
    print("# operator->state overlap Z_i^n (column-normalized |composition|):")
    print("#   rows: 1, O_A, sigma^2, O_22^tt ;  cols: state 0..%d" % (nlev - 1))
    print(np.abs(Zc))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    im = ax.imshow(np.abs(Zc), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(nlev))
    ax.set_xticklabels(["state %d\n$aE$=%.3f" % (n, En[n]) for n in range(nlev)], fontsize=9)
    ax.set_yticks(range(4))
    ax.set_yticklabels(ops, fontsize=11)
    for i in range(4):
        for n in range(nlev):
            v = abs(Zc[i, n])
            ax.text(n, i, "%.2f" % v, ha="center", va="center",
                    color="white" if v < 0.6 else "black", fontsize=10)
    ax.set_title(r"FREE L=%d  operator$\to$state overlap $|Z_i^n|$ (col-norm), $t_0$=%d, $t$=%d"
                 % (dc.L, T0, TMEAS))
    fig.colorbar(im, ax=ax, label="normalized overlap")
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/gevp_overlap_heatmap_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
