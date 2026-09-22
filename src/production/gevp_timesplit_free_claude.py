#!/usr/bin/env python3
# gevp_timesplit_free_claude.py  [FREE -- {1,sigma^2,O_A,O_22} enlarged by a Hankel time-shift (pencil-of-functions)]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 SHIFT=2 NSH=2 python3 gevp_timesplit_free_claude.py
#       ENS=free LREF=2 NVDIR=distill_Nv84 SHIFT=2 NSH=2 python3 gevp_timesplit_free_claude.py
#
# Motivation: {1,sigma^2,O_22} gives vacuum + (1,1,1,1)=2E_1 but NOT the two-meson (2m_sigma): the
# two-meson is HEAVIER than 2E_1 (sub-linear dispersion) and sigma^2 is dominated by the lighter 2E_1,
# so its two-meson piece can't define a clean level.  Enlarge the basis WITHOUT new operators via the
# pencil-of-functions (generalized-pencil / time-shift) trick: from the N-op correlator C(t) build the
# NSH*N Hankel matrix  Chat(t)_{ab} = C(t + (a+b)*SHIFT)  (a,b = 0..NSH-1), then GEVP.  The time shifts
# give independent views that can peel the excited two-meson off the 2E_1 ground.
# Correlators identical to gevp_OA_sig22_free_claude.py.  Ref (GPOF): Aubin-Orginos 1010.0202.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = int(os.environ.get("T0", "2"))
DT_DECOMP = int(os.environ.get("DT_DECOMP", "3"))
SHIFT = int(os.environ.get("SHIFT", "2"))
NSH = int(os.environ.get("NSH", "2"))


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


def build_4x4(twin):
    de.CONTACT = 0.5
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, tw = dc.load_peram(dc.KS[0])
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(tw):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(tw)]
    PA = [Phi[a] @ tt[a, a] for a in range(tw)]
    Q2, clusters, Esh = shell_projector(tau, tw, Nv, DT_DECOMP)
    PQ = [Q2 @ tt[a, a] for a in range(tw)]
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(tw)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(tw)])
    oQ = np.mean([(-np.trace(PQ[a] @ tt[a, a])).real for a in range(tw)])
    CAA = np.zeros(tw)
    C22 = np.zeros(tw)
    CQQ = np.zeros(tw)
    C2A = np.zeros(tw)
    C2Q = np.zeros(tw)
    CAQ = np.zeros(tw)
    for dt in range(tw):
        ns = tw - dt
        aAA = a22 = aQQ = a2A = a2Q = aAQ = 0.0
        for s in range(ns):
            t = s + dt
            aAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
            aQQ += (-np.trace(PQ[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            aAQ += (-np.trace(PA[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, t)).sum().real
            a2A += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            a2Q += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PQ[t] @ tau[t, s])).real
        CAA[dt] = aAA / ns
        C22[dt] = a22 / ns
        CQQ[dt] = aQQ / ns
        C2A[dt] = a2A / ns
        C2Q[dt] = a2Q / ns
        CAQ[dt] = aAQ / ns
    Cts = np.zeros((tw, 4, 4))
    for dt in range(tw):
        Cts[dt] = np.array([
            [1.0,   oA,               o2,               oQ],
            [oA,    CAA[dt] + oA ** 2, C2A[dt] + o2 * oA, CAQ[dt] + oA * oQ],
            [o2,    C2A[dt] + o2 * oA, C22[dt],           C2Q[dt] + o2 * oQ],
            [oQ,    CAQ[dt] + oA * oQ, C2Q[dt] + o2 * oQ, CQQ[dt] + oQ ** 2]])
    return Cts, tw, clusters


def hankel(Cts, shift, nsh):
    twin, N, _ = Cts.shape
    tmax = twin - shift * 2 * (nsh - 1)             # max block offset is (a+b)_max*shift = 2*(nsh-1)*shift
    Big = np.zeros((tmax, nsh * N, nsh * N))
    for t in range(tmax):
        for a in range(nsh):
            for b in range(nsh):
                Big[t, a * N:(a + 1) * N, b * N:(b + 1) * N] = Cts[t + (a + b) * shift]
    return Big


def gevp(Cts, t0, tol=1e-10):
    twin = Cts.shape[0]
    C0 = 0.5 * (Cts[t0] + Cts[t0].T)
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
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s FREE L=%d  {1,sigma^2,O_A,O_22} + time-split (SHIFT=%d NSH=%d)  T0=%d"
          % (tag, dc.L, SHIFT, NSH, T0))
    print("# m_sig~%.3f  2m_sig(two-meson)~%.3f  (2,2)=2E_1~0.52" % (msig, 2 * msig))
    Cts, twin, clusters = build_4x4(twin=None) if False else build_4x4(0)
    print("# shell degeneracies: %s" % [c[1] - c[0] + 1 for c in clusters])
    Big = hankel(Cts, SHIFT, NSH)
    lam, nlev = gevp(Big, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# enlarged levels kept = %d (of %d)" % (nlev, NSH * 4))
    show = min(nlev, 5)
    hdr = "  ".join("m%d" % i for i in range(show))
    print("\n#  t | %s   [look for 2E_1=0.52 and 2m_sig=%.3f plateaus]" % (hdr, 2 * msig))
    for t in range(T0 + 1, Big.shape[0] - 2):
        r = []
        for i in range(show):
            ok = np.isfinite(em[t, i]) and (lam[t, i] * lam[t + 1, i] > 0)
            r.append("%7.4f" % em[t, i] if ok else "  ---  ")
        print("#  %2d | %s" % (t, "  ".join(r)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, Big.shape[0] - 1)
    fig, ax = plt.subplots(figsize=(9, 6))
    cols = ["tab:gray", "tab:blue", "tab:green", "tab:orange", "tab:purple", "tab:brown"]
    for i in range(show):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % len(cols)], marker="o", ms=3, lw=1, label="level %d" % i)
    for y, lab in [(0.0, "vac"), (0.52, r"$(2,2){=}2E_1$"), (2 * msig, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
        ax.text(Big.shape[0] - 3, y + 0.01, lab, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  time-split $\{1,\sigma^2,O_A,O_{22}\}$ GEVP (SHIFT=%d,NSH=%d)"
                 % (dc.L, SHIFT, NSH))
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/gevp_timesplit_free_L%d_s%d_claude.png" % (dc.L, SHIFT)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
