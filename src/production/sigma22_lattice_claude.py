#!/usr/bin/env python3
# sigma22_lattice_claude.py
#   Build the lattice (2,2)=2E1 operator by the PHYSICAL (frame-free) shell projector and confirm 2E1.
#
#   The distillation basis V is Wilson low modes (distill_peram_claude.cu chunk 1); its eigenvalue
#   ordering does NOT match the overlap Dirac shells (at L1 it clusters 8,4,12, not 4,8,12).  But the
#   physical single-fermion shells live in the TIME-dependence of the overlap perambulator:
#     K(dt) = <tau(s+dt, s)>_s = sum_shell Q_shell exp(-E_shell dt) ,  Q_shell = V^dag Pi_shell V,
#   with Q_shell ORTHOGONAL projectors (Q_a Q_b = delta_ab Q_a) since V is t-independent (free config).
#   So eigendecomposing K(dt) recovers the shells: eigenvalue mu_a = exp(-E_a dt), eigenvectors = shells.
#   The (2,2) operator is the l=0 trace over the lambda=2 shell -> vertex Q_{lambda2} (the 2nd-energy block):
#     C_22(dt) = <-Tr[ Q2 tau(s+dt,s) Q2 tau(s,s+dt) ]>_s   -> effmass 2E1.
#   Reference: sigma_00 (full l=0 scalar) -> 2E0.
#
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 sigma22_lattice_claude.py  (or LREF=2 NVDIR=distill_Nv84)

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

DT_DECOMP = int(os.environ.get("DT_DECOMP", "3"))     # dt used to eigendecompose K (shell separation)


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    e0 = 0.5 * msig
    print("# ENS=%s  FREE L=%d  lattice (2,2)=2E1 via physical shell projector  2E0=%.3f  E0~%.3f"
          % (tag, dc.L, msig, e0))

    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = V.shape[1]

    # V t-independence check (free): overlap of V[tsrc0] and V[tsrc0+1] column spaces
    A = V[tsrc0].conj() @ V[tsrc0 + 1].T                # (2Ns,2Ns)? V is (Nt,Nv,2Ns); use mode rows
    vdev = np.abs(np.abs(V[tsrc0]) - np.abs(V[tsrc0 + 1])).max()
    print("# |V[t]|-|V[t+1]| max dev = %.2e (should be ~0 for free, t-independent basis)" % vdev)

    # single-fermion transfer K(dt) = mean_s tau(s+dt, s)
    def Kof(dt):
        ns = twin - dt
        acc = np.zeros((Nv, Nv), complex)
        for s in range(ns):
            acc += tau[s + dt, s]
        return acc / ns

    K = Kof(DT_DECOMP)
    mu, R = np.linalg.eig(K)
    Rinv = np.linalg.inv(R)
    E = -np.log(np.abs(mu)) / DT_DECOMP                 # single-fermion energies
    order = np.argsort(E)
    E = E[order]
    mu = mu[order]
    R = R[:, order]
    Rinv = Rinv[order, :]
    print("# single-fermion shells from K(dt=%d) eigen-decomposition (E = -ln|mu|/dt):" % DT_DECOMP)
    # cluster by E
    clusters = []
    i = 0
    while i < Nv:
        j = i
        while j + 1 < Nv and abs(E[j + 1] - E[i]) < 0.02:
            j += 1
        clusters.append((i, j))
        print("#   shell: E=%.4f  degeneracy=%d  (modes %d..%d)" % (E[i], j - i + 1, i, j))
        i = j + 1
    # lambda=1 = clusters[0] (E~E0), lambda=2 = clusters[1] (E~E1)
    i2, j2 = clusters[1]
    sel = list(range(i2, j2 + 1))
    Q2 = R[:, sel] @ Rinv[sel, :]                       # spectral projector onto the lambda=2 shell
    idem = np.abs(Q2 @ Q2 - Q2).max()
    print("# lambda=2 shell: E1=%.4f (2E1=%.4f), degeneracy=%d, projector idempotency |Q2^2-Q2|=%.2e"
          % (E[i2], 2 * E[i2], len(sel), idem))

    def twopt_vertex(P_snk, P_src):
        C = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                acc += (-np.trace(P_snk @ tau[t, s] @ P_src @ tau[s, t])).real
            C[dt] = acc / ns
        return C

    # sigma_00 uses the physical area*Y00 vertex Phi_00(t); (2,2) uses Q2 (frame-free, per-t identical)
    Phi = [V[tsrc0 + a].T.conj().T @ (w00[:, None] * V[tsrc0 + a].T) for a in range(twin)]
    def twopt_phi(Pl):
        C = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                acc += (-np.trace(Pl[t] @ tau[t, s] @ Pl[s] @ tau[s, t])).real
            C[dt] = acc / ns
        return C
    C00 = twopt_phi(Phi)
    C22 = twopt_vertex(Q2, Q2)

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))
    e22, e00 = eff(C22), eff(C00)
    print("\n  dt |  C_22          m_eff(2,2)  | m_eff(sigma_00)")
    for dt in range(1, twin - 1):
        print("  %2d |  %+.4e     %6.3f     |   %6.3f" % (dt, C22[dt], e22[dt], e00[dt]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(C22, "tab:purple", "o", r"$(2,2)$ shell-projected $=2E_1$"),
                            (C00, "tab:gray", "s", r"$\sigma_{00}$ $=2E_0$")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=4, lw=1.2, label=lab)
    ax.axhline(msig, color="tab:gray", ls="--", lw=1.0, alpha=0.7, label=r"$2E_0=%.3f$" % msig)
    ax.axhline(2 * E[i2], color="tab:purple", ls=":", lw=1.0, alpha=0.7, label=r"$2E_1=%.3f$" % (2 * E[i2]))
    ax.set_ylim(0.2, 1.0)
    ax.set_xlabel("dt")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  lattice $(2,2)$ shell-projected: $2E_1$" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma22_lattice_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
