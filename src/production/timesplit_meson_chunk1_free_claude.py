#!/usr/bin/env python3
# timesplit_meson_chunk1_free_claude.py  [chunk 1: split vertices + split x split correlator, single-op sanity]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 timesplit_meson_chunk1_free_claude.py
#
# Time-split mesonic ops (psibar at t, psi at t+1), kernel = one-step perambulator tau(t,t+1) (contact-free):
#   O_A^s(a)  = psibar(a) [Phi(a) tau(a,a+1)] psi(a+1)     vertex M_A(a) = Phi(a) @ tau(a,a+1)
#   O_22^s(a) = psibar(a) [Q2   tau(a,a+1)] psi(a+1)       vertex M_Q(a) = Q2   @ tau(a,a+1)
# split x split correlator (source (s,s+1), sink (t,t+1), dt=t-s):
#   C^s(dt) = mean_s  -Tr[ M(t) tau(t+1,s+1) M(s)^dag tau(s,t) ]
# Sanity: O_22^s alone should give a clean pair plateau; compare to equal-time O_22.  See
# timesplit_meson_impl_plan_claude.md.  Q2 = project-specific lambda=2 spectral projector (sigma22_lattice).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

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
    R = R[:, order]
    Rinv = Rinv[order, :]
    E = E[order]
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


def effmass(C):
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


def split_corr(M, tau, twin):
    # C^s(dt) = mean_s -Tr[ M[t] tau[t+1,s+1] M[s]^dag tau[s,t] ],  t=s+dt ; needs t+1 <= twin-1
    C = np.full(twin, np.nan)
    Md = [m.conj().T for m in M]
    for dt in range(twin - 1):
        acc = 0.0
        cnt = 0
        for s in range(twin - 1 - dt):
            t = s + dt
            acc += (-np.trace(M[t] @ tau[t + 1, s + 1] @ Md[s] @ tau[s, t])).real
            cnt += 1
        if cnt > 0:
            C[dt] = acc / cnt
    return C


def eqtime_QQ(PQ, tau, twin):
    C = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        acc = 0.0
        for s in range(ns):
            t = s + dt
            acc += (-np.trace(PQ[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
        C[dt] = acc / ns
    return C


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s FREE L=%d  time-split meson CHUNK 1 (split vertices + split x split)" % (tag, dc.L))
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
    print("# shell degeneracies: %s  (Q2: tr=%.3f, idempotency |Q2^2-Q2|=%.2e)"
          % ([c[1] - c[0] + 1 for c in clusters], np.trace(Q2).real, np.abs(Q2 @ Q2 - Q2).max()))

    # split vertices M(a) = K @ tau(a,a+1)
    MA = [Phi[a] @ tau[a, a + 1] for a in range(twin - 1)] + [np.zeros((Nv, Nv), complex)]
    MQ = [Q2 @ tau[a, a + 1] for a in range(twin - 1)] + [np.zeros((Nv, Nv), complex)]
    # equal-time O_22 vertex for comparison
    PQ = [Q2 @ tt[a, a] for a in range(twin)]

    CAs = split_corr(MA, tau, twin)
    CQs = split_corr(MQ, tau, twin)
    CQe = eqtime_QQ(PQ, tau, twin)

    emAs = effmass(CAs)
    emQs = effmass(CQs)
    emQe = effmass(CQe)
    print("\n#  t |   O_A^s        O_22^s     |  m_A^s    m_22^s   m_22(eqtime)")
    for t in range(1, min(twin - 2, 26)):
        def fmt(em, C):
            ok = np.isfinite(em[t]) and np.isfinite(C[t]) and np.isfinite(C[t + 1]) and C[t] * C[t + 1] > 0
            return "%8.4f" % em[t] if ok else "  ----  "
        print("#  %2d | %10.3e  %10.3e |  %s  %s  %s"
              % (t, CAs[t], CQs[t], fmt(emAs, CAs), fmt(emQs, CQs), fmt(emQe, CQe)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(twin - 1)
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    for em, C, col, mk, lab in [(emAs, CAs, "tab:red", "o", r"$O_A^s$ (split)"),
                                 (emQs, CQs, "tab:blue", "s", r"$O_{22}^s$ (split)"),
                                 (emQe, CQe, "tab:cyan", "^", r"$O_{22}$ (eq-time)")]:
        g = np.isfinite(em) & np.isfinite(C[:-1]) & np.isfinite(C[1:]) & (C[:-1] * C[1:] > 0)
        ax.plot(ts[g], em[g], color=col, marker=mk, ms=4, lw=1, label=lab)
    for y, l in [(0.52, r"$(2,2){=}2E_1$"), (0.756, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls=":", lw=0.8, alpha=0.4)
        ax.text(twin - 3, y + 0.008, l, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  time-split single-op sanity (chunk 1)" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/timesplit_meson_chunk1_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
