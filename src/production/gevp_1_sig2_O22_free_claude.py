#!/usr/bin/env python3
# gevp_1_sig2_O22_free_claude.py  [FREE -- minimal {1, sigma^2, O_22} GEVP; drops O_A]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 gevp_1_sig2_O22_free_claude.py
#       ENS=free LREF=2 NVDIR=distill_Nv84 python3 gevp_1_sig2_O22_free_claude.py
#
# Minimal basis for the three wanted states:
#   1     -> vacuum
#   O_22 = psibar Q2 tilde_tau psi  -> genuine (1,1,1,1)=2E_1   (clean lambda=2 projector)
#   sigma^2 = sigma_00^2            -> two-meson 2m_sigma        (disconnected piece)
# O_A is DROPPED: it injects a spurious level that decays to null (its dominant lambda=1 piece is not a
# genuine state), polluting the ground GEVP level (gevp_OA_O22_2op_free_claude.py + shell decomp).
# Vertices/contractions identical to gevp_OA_sig22_free_claude.py.  FULL correlators; identity = vacuum.

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


def gevp(Cts, t0, tol=1e-11):
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


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s FREE L=%d  {1, sigma^2, O_22} GEVP  m_sig~%.3f 2m_sig~%.3f (2,2)~0.52  T0=%d"
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

    Q2, clusters, Esh = shell_projector(tau, twin, Nv, DT_DECOMP)
    print("# shell degeneracies (should be 4,8,12): %s" % [c[1] - c[0] + 1 for c in clusters])
    PQ = [Q2 @ tt[a, a] for a in range(twin)]

    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    oQ = np.mean([(-np.trace(PQ[a] @ tt[a, a])).real for a in range(twin)])
    print("# one-points: <sigma^2>=%.4e  <O_22>=%.4e" % (o2, oQ))

    C22 = np.zeros(twin)
    CQQ = np.zeros(twin)
    C2Q = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        a22 = 0.0
        aQQ = 0.0
        a2Q = 0.0
        for s in range(ns):
            t = s + dt
            aQQ += (-np.trace(PQ[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum().real
            a2Q += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PQ[t] @ tau[t, s])).real
        C22[dt] = a22 / ns
        CQQ[dt] = aQQ / ns
        C2Q[dt] = a2Q / ns

    # FULL 3x3; order {0:1, 1:sigma^2, 2:O_22}
    Cts = np.zeros((twin, 3, 3))
    for dt in range(twin):
        Cts[dt] = np.array([
            [1.0,   o2,               oQ],
            [o2,    C22[dt],          C2Q[dt] + o2 * oQ],
            [oQ,    C2Q[dt] + o2 * oQ, CQQ[dt] + oQ ** 2]])

    lam, vecs, nlev = gevp(Cts, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# levels kept = %d" % nlev)
    print("\n#  t |  m0(vac)   m1        m2     | m1 eigvec[1,s2,O22]           m2 eigvec")
    for t in range(T0 + 1, twin - 2):
        r = []
        for i in range(3):
            ok = i < nlev and np.isfinite(em[t, i]) and (lam[t, i] * lam[t + 1, i] > 0)
            r.append("%8.4f" % em[t, i] if ok else "  ---   ")
        def vs(k):
            if k >= nlev:
                return "        --        "
            v = vecs[t][:, k]
            v = v / (np.linalg.norm(v) + 1e-300)
            return "[%+.2f %+.2f %+.2f]" % (v[0].real, v[1].real, v[2].real)
        print("#  %2d | %s  %s  %s | %s %s" % (t, r[0], r[1], r[2], vs(1), vs(2)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:green", "tab:blue"]
    labs = ["level 0 (vac)", "level 1", "level 2"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nlev):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % 3], marker="o", ms=3, lw=1, label=labs[i])
    for y, lab in [(0.0, "vac"), (msig, r"$m_\sigma$"), (0.52, r"$(2,2){=}2E_1$"), (2 * msig, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
        ax.text(twin - 3, y + 0.01, lab, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  $\{1,\ \sigma_{00}^2,\ O_{22}\}$ GEVP levels" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/gevp_1_sig2_O22_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
