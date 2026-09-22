#!/usr/bin/env python3
# gevp_twosigma_split_free_claude.py  [time-split two-sigma operator; {1,sigma^2,O_2sigma} then +O_22]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 SPLIT=1 python3 gevp_twosigma_split_free_claude.py
#       ENS=free LREF=2 NVDIR=distill_Nv84 SPLIT=1 python3 gevp_twosigma_split_free_claude.py
#
# Two-meson lives in the FOUR-fermion sector; a bilinear cannot reach it (free particle number).  The
# genuine time-split two-meson interpolator is  O_2sigma(t) = sigma_00(t) sigma_00(t+SPLIT), a second
# four-fermion op independent of sigma^2 = sigma_00(t)^2.
# Single-sigma meson propagator (decays as m_sigma = 2E0):  M(t,s) = -Tr[Phi(t) tau(t,s) Phi(s) tau(s,t)].
# Two-meson (factorized, the CS^2 piece with split times):
#   <O_2s(t) O_2s(0)>_c = mean_s [ M(t,s) M(t+d,s+d) + M(t,s+d) M(t+d,s) ]
#   <sigma^2(t) O_2s(0)>_c = mean_s  2 M(t,s) M(t+d,s)      (sigma^2 both at s)
# O_22 (bilinear, 2E1) is added in the +O_22 run; its cross with O_2s is 4f-vs-2f connected (small) -> set 0.
# See timesplit_meson_impl_plan_claude.md.

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
SPLIT = int(os.environ.get("SPLIT", "1"))
ADD22 = int(os.environ.get("ADD22", "1"))
BASIS22S = int(os.environ.get("BASIS22S", "0"))    # 1 -> clean {1, O_22, O_2sigma} (drop sigma^2)
BASIS5 = int(os.environ.get("BASIS5", "0"))        # 1 -> {1, sigma^2, O_2sigma, O_A, O_22}


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
    i2, j2 = None, None
    Ecl = E[order]
    clusters = []
    i = 0
    while i < Nv:
        j = i
        while j + 1 < Nv and abs(Ecl[j + 1] - Ecl[i]) < 0.02:
            j += 1
        clusters.append((i, j))
        i = j + 1
    i2, j2 = clusters[1]
    sel = list(range(i2, j2 + 1))
    Q2 = R[:, sel] @ Rinv[sel, :]
    return Q2, clusters


def effmass(C):
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


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
    d = SPLIT
    print("# ENS=%s FREE L=%d  two-sigma split d=%d  {1,sigma^2,O_2sigma%s}  T0=%d"
          % (tag, dc.L, d, ",O_22" if ADD22 else "", T0))
    print("# m_sig~%.3f  2m_sig(two-meson)~%.3f  (2,2)=2E_1~0.52" % (msig, 2 * msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    # single-sigma meson propagator M(t,s) = -Tr[Phi(t) tau(t,s) Phi(s) tau(s,t)]
    Msig = np.zeros((twin, twin))
    for t in range(twin):
        for s in range(twin):
            Msig[t, s] = (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
    emM = effmass(np.array([Msig[k, 0] for k in range(twin)]))
    print("# single-sigma M(t,0) effmass (expect ~m_sig=%.3f):  %s"
          % (msig, "  ".join("%.3f" % emM[k] for k in range(3, min(12, twin - 1)))))

    osig = np.mean([(-np.trace(Phi[a] @ tt[a, a])).real for a in range(twin)])
    o2 = np.mean([(-np.trace(Phi[a] @ tt[a, a] @ Phi[a] @ tt[a, a])).real for a in range(twin)])
    o2s = np.mean([osig ** 2 + Msig[a, a + d] for a in range(twin - d)])

    hi = twin - 1
    C22 = np.full(twin, np.nan)
    C2s = np.full(twin, np.nan)          # <sigma^2 O_2s>_c
    Css = np.full(twin, np.nan)          # <O_2s  O_2s>_c
    CQQ = np.full(twin, np.nan)
    C2Q = np.full(twin, np.nan)
    CAA = np.full(twin, np.nan)          # O_A x O_A
    CAQ = np.full(twin, np.nan)          # O_A x O_22
    C2A = np.full(twin, np.nan)          # sigma^2 x O_A
    Q2, clusters = shell_projector(tau, twin, Nv, DT_DECOMP)
    PQ = [Q2 @ tt[a, a] for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    for dt in range(twin):
        ns = twin - dt
        a22 = 0.0
        for s in range(ns):
            a22 += (dc.W10 * de.diags_pair(Phi, tau, s, s + dt)).sum().real
        C22[dt] = a22 / ns
        # two-sigma pieces (need t+d, s+d <= hi)
        acc2s = 0.0
        accss = 0.0
        accQQ = 0.0
        acc2Q = 0.0
        accAA = 0.0
        accAQ = 0.0
        accA2 = 0.0
        cnt = 0
        for s in range(ns):
            t = s + dt
            if t + d > hi or s + d > hi:
                continue
            acc2s += 2.0 * Msig[t, s] * Msig[t + d, s]
            accss += Msig[t, s] * Msig[t + d, s + d] + Msig[t, s + d] * Msig[t + d, s]
            accQQ += (-np.trace(PQ[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            acc2Q += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PQ[t] @ tau[t, s])).real
            accAA += (-np.trace(PA[t] @ tau[t, s] @ PA[s] @ tau[s, t])).real
            accAQ += (-np.trace(PA[t] @ tau[t, s] @ PQ[s] @ tau[s, t])).real
            accA2 += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ PA[t] @ tau[t, s])).real
            cnt += 1
        if cnt > 0:
            C2s[dt] = acc2s / cnt
            Css[dt] = accss / cnt
            CQQ[dt] = accQQ / cnt
            C2Q[dt] = acc2Q / cnt
            CAA[dt] = accAA / cnt
            CAQ[dt] = accAQ / cnt
            C2A[dt] = accA2 / cnt

    emS = effmass(Css)
    print("# O_2sigma <O_2s O_2s>_c effmass (expect ~2m_sig=%.3f):  %s"
          % (2 * msig, "  ".join("%.3f" % emS[k] if np.isfinite(emS[k]) else "--" for k in range(3, min(14, twin - 1)))))

    oQ = np.mean([(-np.trace(PQ[a] @ tt[a, a])).real for a in range(twin)])
    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    if BASIS5:
        # {0:1, 1:sigma^2, 2:O_2sigma, 3:O_A, 4:O_22}
        N = 5
        Cts = np.full((twin, N, N), np.nan)
        for dt in range(twin):
            M = np.zeros((N, N))
            M[0, 0] = 1.0
            ops1 = [o2, o2s, oA, oQ]
            for k in range(4):
                M[0, k + 1] = M[k + 1, 0] = ops1[k]
            M[1, 1] = C22[dt]
            M[1, 2] = M[2, 1] = C2s[dt] + o2 * o2s
            M[1, 3] = M[3, 1] = C2A[dt] + o2 * oA
            M[1, 4] = M[4, 1] = C2Q[dt] + o2 * oQ
            M[2, 2] = Css[dt] + o2s ** 2
            M[2, 3] = M[3, 2] = 0.0 + o2s * oA       # O_2s x O_A (4f-2f) set 0
            M[2, 4] = M[4, 2] = 0.0 + o2s * oQ       # O_2s x O_22 (4f-2f) set 0
            M[3, 3] = CAA[dt] + oA ** 2
            M[3, 4] = M[4, 3] = CAQ[dt] + oA * oQ
            M[4, 4] = CQQ[dt] + oQ ** 2
            Cts[dt] = M
    elif BASIS22S:
        # clean CONNECTED 2x2 {0:O_22, 1:O_2sigma} -- no identity (O_2sigma is vacuum-dominated, so the
        # identity-subtracted GEVP is ill-conditioned; use connected correlators directly).  Cross ~0
        # (4f-2f, particle number) so the GEVP is near-diagonal: O_22->2E1, O_2sigma->two-meson.
        N = 2
        Cts = np.full((twin, N, N), np.nan)
        for dt in range(twin):
            M = np.zeros((N, N))
            M[0, 0] = CQQ[dt]
            M[1, 1] = Css[dt]
            M[0, 1] = M[1, 0] = 0.0                  # O_22 x O_2s connected cross set to 0
            Cts[dt] = M
    else:
        # {0:1, 1:sigma^2, 2:O_2sigma [,3:O_22]}
        N = 4 if ADD22 else 3
        Cts = np.full((twin, N, N), np.nan)
        for dt in range(twin):
            M = np.zeros((N, N))
            M[0, 0] = 1.0
            M[0, 1] = M[1, 0] = o2
            M[0, 2] = M[2, 0] = o2s
            M[1, 1] = C22[dt]
            M[1, 2] = M[2, 1] = C2s[dt] + o2 * o2s
            M[2, 2] = Css[dt] + o2s ** 2
            if ADD22:
                M[0, 3] = M[3, 0] = oQ
                M[1, 3] = M[3, 1] = C2Q[dt] + o2 * oQ
                M[2, 3] = M[3, 2] = 0.0 + o2s * oQ      # O_2s x O_22 (4f-2f connected) set to 0
                M[3, 3] = CQQ[dt] + oQ ** 2
            Cts[dt] = M

    good = np.array([np.all(np.isfinite(Cts[t])) for t in range(twin)])
    tvec = np.where(good)[0]
    Cg = Cts[tvec]
    t0i = int(np.argmin(np.abs(tvec - T0)))
    lam, nlev = gevp(Cg, t0i)
    em = effmass(lam)
    print("# levels kept = %d of %d" % (nlev, N))
    print("\n#  t | " + "  ".join("m%d" % i for i in range(nlev)) + "   [vac, 2E_1=0.52, 2m_sig=%.3f]" % (2 * msig))
    for i in range(len(tvec) - 2):
        t = tvec[i]
        if t <= T0:
            continue
        r = []
        for k in range(nlev):
            ok = np.isfinite(em[i, k]) and (lam[i, k] * lam[i + 1, k] > 0)
            r.append("%7.4f" % em[i, k] if ok else "  ---  ")
        print("#  %2d | %s" % (t, "  ".join(r)))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = tvec[:-1]
    fig, ax = plt.subplots(figsize=(9, 6))
    cols = ["tab:gray", "tab:green", "tab:red", "tab:blue"]
    for k in range(nlev):
        g = np.isfinite(em[:, k]) & (lam[:-1, k] * lam[1:, k] > 0)
        ax.plot(ts[g], em[g, k], color=cols[k % 4], marker="o", ms=3, lw=1, label="level %d" % k)
    for y, l in [(0.0, "vac"), (0.52, r"$(2,2){=}2E_1$"), (2 * msig, r"$2m_\sigma$")]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
        ax.text(twin - 3, y + 0.01, l, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  two-$\sigma$ split (d=%d) $\{1,\sigma^2,O_{2\sigma}%s\}$"
                 % (dc.L, d, r",O_{22}" if ADD22 else ""))
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    suff = "_5op" if BASIS5 else ("_22s" if BASIS22S else ("_q" if ADD22 else ""))
    out = "figs/gevp_twosigma_split_free_L%d_d%d%s_claude.png" % (dc.L, d, suff)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
