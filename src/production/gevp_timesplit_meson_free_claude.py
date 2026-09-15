#!/usr/bin/env python3
# gevp_timesplit_meson_free_claude.py  [chunk 2: full 6x6 {1,sigma^2,O_A,O_22,O_A^s,O_22^s} GEVP]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 gevp_timesplit_meson_free_claude.py
#       ENS=free LREF=2 NVDIR=distill_Nv84 python3 gevp_timesplit_meson_free_claude.py
#
# Equal-time ops (vertex at time a): O_A = psibar Phi tilde_tau psi ; O_22 = psibar Q2 tilde_tau psi.
# Time-split ops (psibar at a, psi at a+1, kernel = tau(a,a+1), contact-free):
#   O_A^s = psibar(a) [Phi(a) tau(a,a+1)] psi(a+1) ; O_22^s = psibar(a) [Q2 tau(a,a+1)] psi(a+1).
# Bilinear cross rule (sink O_x at t with psi-leg t+dx, source O_y at s with psi-leg s+dy):
#   C_xy(dt) = mean_s  -Tr[ Kx(t) tau(t+dx, s+dy) Ky(s)^dag tau(s,t) ]    (dx,dy in {0 eqtime, 1 split})
# sigma^2 x bilinear cross (sigma^2 at source s, bilinear at sink t, psi-leg t+dx):
#   C_2x(dt) = mean_s  -2 Tr[ Phi(s) tilde_tau(s) Phi(s) tau(s,t) Kx(t) tau(t+dx, s) ]
# one-point <O_x> = mean_a -Tr[ Kx(a) tau(a+dx, a) ].   See timesplit_meson_impl_plan_claude.md.

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


def bilinear_cross(Kx, dx, Ky, dy, tau, twin):
    # C_xy(dt) = mean_s -Tr[ Kx[t] tau[t+dx, s+dy] Ky[s]^dag tau[s,t] ] ; needs t+dx, s+dy <= twin-1
    C = np.full(twin, np.nan)
    Kyd = [k.conj().T for k in Ky]
    hi = twin - 1
    for dt in range(twin):
        acc = 0.0
        cnt = 0
        smax = twin - dt
        for s in range(smax):
            t = s + dt
            if t + dx > hi or s + dy > hi:
                continue
            acc += (-np.trace(Kx[t] @ tau[t + dx, s + dy] @ Kyd[s] @ tau[s, t])).real
            cnt += 1
        if cnt > 0:
            C[dt] = acc / cnt
    return C


def sig2_cross(Phi, tt, Kx, dx, tau, twin):
    # C_2x(dt) = mean_s -2 Tr[ Phi[s] tt[s,s] Phi[s] tau[s,t] Kx[t] tau[t+dx, s] ]
    C = np.full(twin, np.nan)
    hi = twin - 1
    for dt in range(twin):
        acc = 0.0
        cnt = 0
        for s in range(twin - dt):
            t = s + dt
            if t + dx > hi:
                continue
            acc += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Kx[t] @ tau[t + dx, s])).real
            cnt += 1
        if cnt > 0:
            C[dt] = acc / cnt
    return C


def onepoint(Kx, dx, tau, twin):
    hi = twin - 1
    vals = []
    for a in range(twin):
        if a + dx > hi:
            continue
        vals.append((-np.trace(Kx[a] @ tau[a + dx, a])).real)
    return float(np.mean(vals))


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
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s FREE L=%d  {1,sigma^2,O_A,O_22,O_A^s,O_22^s} 6x6 GEVP  T0=%d" % (tag, dc.L, T0))
    print("# 2m_sig(two-meson)~%.3f  (2,2)=2E_1~0.52" % (2 * msig))
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
    print("# shell degeneracies: %s" % [c[1] - c[0] + 1 for c in clusters])

    # bilinear vertices: order A, 22, A^s, 22^s  with psi-leg shift d
    Zc = np.zeros((Nv, Nv), complex)
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    PQ = [Q2 @ tt[a, a] for a in range(twin)]
    MA = [Phi[a] @ tau[a, a + 1] if a < twin - 1 else Zc for a in range(twin)]
    MQ = [Q2 @ tau[a, a + 1] if a < twin - 1 else Zc for a in range(twin)]
    bil = [("O_A", PA, 0), ("O_22", PQ, 0), ("O_A^s", MA, 1), ("O_22^s", MQ, 1)]

    # one-points (identity row); order in GEVP: {0:1, 1:sigma^2, 2:O_A, 3:O_22, 4:O_A^s, 5:O_22^s}
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    ob = [onepoint(K, d, tau, twin) for (_, K, d) in bil]
    print("# one-points: sig2=%.3e  O_A=%.3e O_22=%.3e O_A^s=%.3e O_22^s=%.3e"
          % (o2, ob[0], ob[1], ob[2], ob[3]))

    # bilinear 4x4 block (connected)
    Cbb = {}
    for a in range(4):
        for b in range(a, 4):
            C = bilinear_cross(bil[a][1], bil[a][2], bil[b][1], bil[b][2], tau, twin)
            Cbb[(a, b)] = C
            Cbb[(b, a)] = C
    # sigma^2 diagonal + sigma^2 x bilinear
    C22 = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        acc = 0.0
        for s in range(ns):
            acc += (dc.W10 * de.diags_pair(Phi, tau, s, s + dt)).sum().real
        C22[dt] = acc / ns
    C2b = [sig2_cross(Phi, tt, bil[a][1], bil[a][2], tau, twin) for a in range(4)]

    # assemble 6x6 FULL correlator (identity carries vacuum)
    obf = [o2] + ob                                   # one-points for {sigma^2, O_A, O_22, O_A^s, O_22^s}
    N = 6
    Cts = np.full((twin, N, N), np.nan)
    for dt in range(twin):
        M = np.zeros((N, N))
        M[0, 0] = 1.0
        for k in range(5):
            M[0, k + 1] = obf[k]
            M[k + 1, 0] = obf[k]
        M[1, 1] = C22[dt]
        for a in range(4):
            v = C2b[a][dt]
            M[1, a + 2] = v + o2 * ob[a]
            M[a + 2, 1] = v + o2 * ob[a]
        for a in range(4):
            for b in range(4):
                M[a + 2, b + 2] = Cbb[(a, b)][dt] + ob[a] * ob[b]
        Cts[dt] = M

    # drop times with any nan (edge from split shifts)
    good_t = np.array([np.all(np.isfinite(Cts[t])) for t in range(twin)])
    tvec = np.where(good_t)[0]
    Cg = Cts[tvec]
    # map T0 to index within good times
    t0i = int(np.argmin(np.abs(tvec - T0)))
    lam, nlev = gevp(Cg, t0i)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
    print("# levels kept = %d (of 6)" % nlev)
    show = min(nlev, 6)
    print("\n#  t | " + "  ".join("m%d" % i for i in range(show)) + "   [vac, 2E_1=0.52, 2m_sig=%.3f]" % (2 * msig))
    for i in range(len(tvec) - 2):
        t = tvec[i]
        if t <= T0:
            continue
        r = []
        for k in range(show):
            ok = np.isfinite(em[i, k]) and (lam[i, k] * lam[i + 1, k] > 0)
            r.append("%7.4f" % em[i, k] if ok else "  ---  ")
        print("#  %2d | %s" % (t, "  ".join(r)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = tvec[:-1]
    fig, ax = plt.subplots(figsize=(9, 6))
    cols = ["tab:gray", "tab:blue", "tab:green", "tab:orange", "tab:purple", "tab:brown"]
    for k in range(show):
        g = np.isfinite(em[:, k]) & (lam[:-1, k] * lam[1:, k] > 0)
        ax.plot(ts[g], em[g, k], color=cols[k % len(cols)], marker="o", ms=3, lw=1, label="level %d" % k)
    for y, l in [(0.0, "vac"), (0.52, r"$(2,2){=}2E_1$"), (2 * msig, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
        ax.text(twin - 3, y + 0.01, l, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  time-split 6x6 $\{1,\sigma^2,O_A,O_{22},O_A^s,O_{22}^s\}$ GEVP" % dc.L)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/gevp_timesplit_meson_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
