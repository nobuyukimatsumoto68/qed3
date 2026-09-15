#!/usr/bin/env python3
# ell2_gevp_free_claude.py
#   Sharpen the E0+E1 plateau with a 2-operator l=2 (icosahedral H) GEVP.
#     O1_2M = psibar Y_2M psi              (local density, vertex Phi_2M)
#     O2_2M = psibar Y_2M tilde_tau psi    (smeared/bilocal, vertex Phi_2M @ tilde_tau) -- analog of O_A
#   Correlator matrix, M-summed (rotation-invariant l=2 channel):
#     C_ij(dt) = sum_M mean_s [ -Tr( Pi_2M(t) tau(t,s) Pj_2M(s) tau(s,t) ) ] ,  t=s+dt.
#   GEVP C(dt) v = lam C(t0) v  (rebased: whiten by the positive part of C(t0)); level 0 -> E0+E1.
#
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 ell2_gevp_free_claude.py   (or LREF=2 NVDIR=distill_Nv84)
#   env: T0 (rebase time, default 4).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import math
import numpy as np
import distill_contract_claude as dc

T0 = int(os.environ.get("T0", "4"))


def real_y2m(xhat):
    x = xhat[:, 0]
    y = xhat[:, 1]
    z = xhat[:, 2]
    c1 = 0.25 * math.sqrt(5.0 / math.pi)
    c2 = 0.5 * math.sqrt(15.0 / math.pi)
    Y = np.zeros((xhat.shape[0], 5))
    Y[:, 0] = c2 * x * y
    Y[:, 1] = c2 * y * z
    Y[:, 2] = c1 * (3.0 * z * z - 1.0)
    Y[:, 3] = c2 * x * z
    Y[:, 4] = 0.5 * c2 * (x * x - y * y)
    return Y


def block(Pi_M, Pj_M, tau, twin):
    # M-summed connected two-point between vertex families Pi_M, Pj_M (each a list over M of (twin,Nv,Nv)).
    C = np.zeros(twin)
    for dt in range(twin):
        ns = twin - dt
        acc = 0.0
        for M in range(5):
            Pi = Pi_M[M]
            Pj = Pj_M[M]
            for s in range(ns):
                t = s + dt
                acc += (-np.trace(Pi[t] @ tau[t, s] @ Pj[s] @ tau[s, t])).real
        C[dt] = acc / ns
    return C


def solve_gevp(c, t0):
    # rebased GEVP: whiten by positive part of symmetrized C(t0), return sorted-desc eigenvalues per t.
    twin = c.shape[0]
    C0 = 0.5 * (c[t0] + c[t0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > 1e-11 * abs(wv).max()
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    ng = int(keep.sum())
    lam = np.full((twin, ng), np.nan)
    for t in range(twin):
        Msym = Uk.T @ (0.5 * (c[t] + c[t].T)) @ Uk
        try:
            lam[t] = np.sort(np.linalg.eigvals(Msym).real)[::-1]
        except Exception:
            pass
    return lam, ng, wv


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    e0e1 = 0.5 * msig + 0.26 if dc.L == 1 else np.nan
    print("# ENS=%s  FREE L=%d  l=2 GEVP {O_A2=Y_2M tt, O_B2=Y_2M tt^2}  2E0=%.3f  T0=%d" % (tag, dc.L, msig, T0))

    dual = dc.dual_areas_from_mesh()
    sites = dc.load_vec3(dc.GEOM + "pts_n%d.dat" % dc.L)
    sites = sites / np.linalg.norm(sites, axis=1, keepdims=True)
    Y2 = real_y2m(sites)
    wl2 = [np.repeat(dual * Y2[:, M], dc.NS) for M in range(5)]

    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Vt = [V[tsrc0 + a].T for a in range(twin)]
    Phi = [[Vt[a].conj().T @ (wl2[M][:, None] * Vt[a]) for a in range(twin)] for M in range(5)]  # local
    # local (Phi, sigma3-even) is orthogonal to any tilde_tau-smeared op by sigma3-hermiticity, so
    # the two operators are both tilde_tau-family (sigma3-mixed -> they mix): O_A2 = Phi tt, O_B2 = Phi tt^2.
    SmA = [[Phi[M][a] @ tt[a, a] for a in range(twin)] for M in range(5)]                         # O_A2
    SmB = [[SmA[M][a] @ tt[a, a] for a in range(twin)] for M in range(5)]                          # O_B2

    C11 = block(SmA, SmA, tau, twin)
    C22 = block(SmB, SmB, tau, twin)
    C12 = block(SmA, SmB, tau, twin)

    cmat = np.zeros((twin, 2, 2))
    for t in range(twin):
        cmat[t] = np.array([[C11[t], C12[t]], [C12[t], C22[t]]])
    lam, ng, wv = solve_gevp(cmat, T0)
    with np.errstate(all="ignore"):
        emg = np.log(lam[:-1] / lam[1:])
    emg = np.pad(emg, ((0, 0), (0, 2 - ng)), constant_values=np.nan)

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))
    e1, e2 = eff(C11), eff(C22)
    cos0 = C12[T0] / math.sqrt(abs(C11[T0] * C22[T0]))
    print("# C(t0) eigenvalues = %s (kept %d)   cos(O1,O2)@t0 = %.4f"
          % (np.array2string(wv, precision=3), ng, cos0))
    print("\n  dt | m(O_A2)      m(O_B2)     | GEVP lvl0   lvl1")
    for dt in range(1, twin - 2):
        print("  %2d |   %6.3f       %6.3f    |  %6.3f    %6.3f"
              % (dt, e1[dt], e2[dt], emg[dt, 0], emg[dt, 1] if ng > 1 else np.nan))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(C11, "tab:gray", "o", r"$O_{A2}=Y_{2M}\tilde\tau$"),
                            (C22, "tab:olive", "^", r"$O_{B2}=Y_{2M}\tilde\tau^2$")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=3, lw=1, alpha=0.7, label=lab)
    g = np.isfinite(emg[:, 0]) & (lam[:-1, 0] * lam[1:, 0] > 0)
    ax.plot(dts[:len(emg)][g], emg[g, 0], color="tab:red", marker="s", ms=4, lw=1.5,
            label=r"$\ell=2$ GEVP level 0 ($E_0+E_1$)")
    if ng > 1:
        g1 = np.isfinite(emg[:, 1]) & (lam[:-1, 1] * lam[1:, 1] > 0)
        ax.plot(dts[:len(emg)][g1], emg[g1, 1], color="tab:blue", marker="v", ms=4, lw=1.0,
                alpha=0.7, label=r"GEVP level 1")
    ax.axhline(msig, color="tab:gray", ls="--", lw=1.0, alpha=0.6, label=r"$2E_0=%.3f$" % msig)
    if np.isfinite(e0e1):
        ax.axhline(e0e1, color="tab:red", ls=":", lw=1.0, alpha=0.7, label=r"$E_0+E_1\approx%.3f$" % e0e1)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlabel("dt")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  $\ell=2$ 2-op GEVP: $E_0+E_1$ plateau" % dc.L)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/ell2_gevp_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
