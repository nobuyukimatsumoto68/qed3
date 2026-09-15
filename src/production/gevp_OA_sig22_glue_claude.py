#!/usr/bin/env python3
# gevp_OA_sig22_glue_claude.py
#   {1, O_A, sigma^2, O_22^tt} PLAIN GEVP (no rebasing), FIXED t0, in the glueball inv_sqrt_sym convention
#   (glue_gevp_analysis_claude.cu use_rebase=0 branch; inv_sqrt_sym metric from :36):
#     - FIXED metric t0 (whitened ONCE): M(t) = C(t0)^{-1/2} C(t) C(t0)^{-1/2} for all t; dt = t - t0 varies
#     - effmass:   m_eff(t) = log(lambda_s(t)/lambda_s(t+1)) / at   (all ng eigenvalues kept, NO rebase)
#     - overlaps:  weight v = C(t0)^{-1/2} w (unit); overlap Z = C(t0) v (unit); both sign-fixed (dom>0)
#   Rebasing (Vre rotation + top-nstates truncation) is DROPPED per request (truncation subspace is
#   operator-dependent).  Constant (identity) kept in the basis -> vacuum is the m_eff~0 mode.  A FIXED t0
#   metric stays well-conditioned at large t (unlike the moving metric), which matters near the noise floor.
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 gevp_OA_sig22_glue_claude.py
#   env:  T0 (fixed metric time, default 3), TMEAS (overlap readout t, default 12), RTOL (1e-8)

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de
import sigma22_spininsert_test_claude as st

T0 = int(os.environ.get("T0", "3"))                  # FIXED metric time t0 (dt = t - t0 varies)
TMEAS = int(os.environ.get("TMEAS", "12"))           # overlap-readout measurement time t
RTOL = float(os.environ.get("RTOL", "1e-8"))
AT = 1.0                                              # report a_t*m (a_t folded into the lattice units)


def inv_sqrt_sym(M, rtol):
    # regularized C^{-1/2} via symmetric eigendecomposition (glue_gevp_analysis_claude.cu:36)
    d, Vv = np.linalg.eigh(0.5 * (M + M.T))
    dmax = d.max()
    idd = np.where((d > rtol * dmax) & (d > 0.0), 1.0 / np.sqrt(np.abs(d)), 0.0)
    return (Vv * idd) @ Vv.T


def gevp_levels(Csym, t0, at, rtol):
    # PLAIN GEVP (no rebase), FIXED t0: metric C(t0) whitened ONCE, M(t) = C(t0)^{-1/2} C(t) C(t0)^{-1/2}
    # for all t; lambda_s(t) = e^{-E_s (t - t0)}.  Effmass = log(lambda_s(t)/lambda_s(t+1))/at (dt = t - t0
    # varies; the metric is FIXED at t0, so it stays well-conditioned at large t unlike the moving metric).
    si0 = inv_sqrt_sym(Csym[t0], rtol)
    twin = len(Csym)
    ng = Csym[0].shape[0]
    lam = np.full((twin, ng), np.nan)
    for t in range(twin):
        M = si0 @ Csym[t] @ si0
        lam[t] = np.linalg.eigvalsh(0.5 * (M + M.T))  # ascending lambda
    eff = np.full((twin, ng), np.nan)
    for t in range(twin - 1):
        good = (lam[t] > 0) & (lam[t + 1] > 0)
        eff[t, good] = np.log(lam[t, good] / lam[t + 1, good]) / at
    return eff                                        # eff[t, s]: s ascending lambda = descending energy


def gevp_overlaps(Csym, tmeas, t0, rtol):
    # FIXED t0 GEVP at measurement time tmeas: M = C(t0)^{-1/2} C(tmeas) C(t0)^{-1/2}; weight v = C(t0)^{-1/2} w
    # (unit), overlap Z = C(t0) v (unit); state energy E = -log(lambda)/(tmeas - t0).
    si0 = inv_sqrt_sym(Csym[t0], rtol)
    M = si0 @ Csym[tmeas] @ si0
    ev, W = np.linalg.eigh(0.5 * (M + M.T))           # ascending
    C0 = Csym[t0]
    ng = C0.shape[0]
    si = si0                                          # weight uses the fixed-t0 metric
    order = np.argsort(ev)[::-1]                       # descending lambda = ascending energy
    vlist = []
    Zlist = []
    Elist = []
    for col in order:
        w = W[:, col]
        v = si @ w
        v = v / np.linalg.norm(v)
        imax = np.argmax(np.abs(v))
        if v[imax] < 0:
            v = -v
        Z = C0 @ v
        Z = Z / np.linalg.norm(Z)
        if Z[imax] < 0:
            Z = -Z
        vlist.append(v)
        Zlist.append(Z)
        Elist.append(-np.log(ev[col]) / ((tmeas - t0) * AT) if ev[col] > 0 else np.nan)
    return np.array(vlist).T, np.array(Zlist).T, np.array(Elist)   # (ng, ng): columns = states, energy asc


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
        M = np.array([
            [1.0,   oA,               o2,               oQ],
            [oA,    CAA[dt] + oA ** 2, C2A[dt] + o2 * oA, CAQ[dt] + oA * oQ],
            [o2,    C2A[dt] + o2 * oA, C22[dt],           C2Q[dt] + o2 * oQ],
            [oQ,    CAQ[dt] + oA * oQ, C2Q[dt] + o2 * oQ, CQQ[dt] + oQ ** 2]])
        return 0.5 * (M + M.T)
    Csym = [Cmat(dt) for dt in range(twin)]
    ops = ["1", "O_A", r"$\sigma^2$", r"$O_{22}^{\tilde\tau}$"]

    eff = gevp_levels(Csym, T0, AT, RTOL)
    print("# PLAIN GEVP (no rebase), FIXED t0=%d  L=%d  rtol=%.0e  2E0=%.3f 2m=%.3f"
          % (T0, dc.L, RTOL, msig, 2 * msig))
    print("#  t |  E(vac..2mes) ascending   [vac~0; (2,2)~0.52; 2m=%.3f]" % (2 * msig))
    for t in range(T0 + 1, twin - 2):
        vals = np.sort(eff[t][np.isfinite(eff[t])])
        r = "  ".join("%7.4f" % x for x in vals)
        print("  %2d | %s" % (t, r))

    v, Z, En = gevp_overlaps(Csym, TMEAS, T0, RTOL)
    print("\n# overlaps at t=%d (glue convention, Z=C(t0)v col-norm), state energies asc = %s"
          % (TMEAS, np.array2string(En, precision=3)))
    print("# rows: 1, O_A, sigma^2, O_22^tt ; cols: state 0..3")
    print(np.array2string(np.abs(Z), precision=3, suppress_small=True))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # -- levels plot --
    ts = np.arange(twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    cols = ["tab:gray", "tab:red", "tab:green", "tab:blue"]
    effsort = np.sort(eff, axis=1)                     # ascending energy per t (nan last)
    ng = eff.shape[1]
    for s in range(ng):
        y = effsort[:, s]
        g = np.isfinite(y)
        ax.plot(ts[g], y[g], color=cols[s % 4], marker="o", ms=3, lw=1, label="level %d" % s)
    for yv, lab in [(0.0, "vac"), (msig, r"$m_\sigma$"), (0.52, r"$(2,2)$"), (2 * msig, r"$2m_\sigma$")]:
        ax.axhline(yv, color="k", ls="--", lw=0.8, alpha=0.4)
        ax.text(twin - 3, yv + 0.01, lab, fontsize=8, alpha=0.7)
    ax.set_ylim(-0.1, 1.2)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$a_t\,m_\mathrm{eff}(t) = \log[\lambda(t)/\lambda(t{+}1)]$")
    ax.set_title(r"FREE L=%d PLAIN GEVP $\{1,O_A,\sigma^2,O_{22}^{\tilde\tau}\}$ (fixed $t_0$=%d)"
                 % (dc.L, T0))
    ax.legend(fontsize=9, ncol=2)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out1 = "figs/gevp_OA_sig22_glue_free_L%d_claude.png" % dc.L
    fig.savefig(out1, dpi=130)
    plt.close(fig)
    # -- overlap heatmap --
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    im = ax.imshow(np.abs(Z), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["state %d\n$aE$=%.3f" % (n, En[n]) for n in range(4)], fontsize=9)
    ax.set_yticks(range(4))
    ax.set_yticklabels(ops, fontsize=11)
    for i in range(4):
        for n in range(4):
            val = abs(Z[i, n])
            ax.text(n, i, "%.2f" % val, ha="center", va="center",
                    color="white" if val < 0.6 else "black", fontsize=10)
    ax.set_title(r"FREE L=%d plain-GEVP overlap $|Z_i^n|=|C(t_0)v|$ (fixed $t_0$=%d, t=%d)" % (dc.L, T0, TMEAS))
    fig.colorbar(im, ax=ax, label="normalized overlap")
    fig.tight_layout()
    out2 = "figs/gevp_OA_sig22_glue_overlap_L%d_claude.png" % dc.L
    fig.savefig(out2, dpi=130)
    plt.close(fig)
    print("\n# -> %s\n# -> %s" % (out1, out2))


if __name__ == "__main__":
    main()
