#!/usr/bin/env python3
# o_b_test_free_claude.py  -- second (2,2) interpolator O_B = psibar tilde_tau^2 psi.
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 o_b_test_free_claude.py   (or LREF=2 NVDIR=distill_Nv84)
# O_A vertex = Phi_00 @ tilde_tau ; O_B vertex = Phi_00 @ tilde_tau^2 = O_A @ tilde_tau (extra smearing).
# Tests: (a) <O_B O_B> effmass -> does it overlap (2,2)?  (b) <O_A O_B> -> independent?  (c) {O_A,O_B}
# 2x2 GEVP (single-meson sector, t0=4) -> does it sharpen the (2,2) plateau?

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

T0 = int(os.environ.get("T0", "4"))


def main():
    tag = dc.ENS.split("nu0")[0]
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# ENS=%s  FREE L=%d  O_B = psibar tilde_tau^2 psi  m_sig~%.3f  (2,2)=2E1" % (tag, dc.L, msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]              # O_A vertex
    PB = [PA[a] @ tt[a, a] for a in range(twin)]               # O_B vertex = Phi_00 tilde_tau^2

    oA = np.mean([(-np.trace(PA[a] @ tt[a, a])).real for a in range(twin)])
    oB = np.mean([(-np.trace(PB[a] @ tt[a, a])).real for a in range(twin)])
    print("# one-points: <O_A> = %.4e   <O_B> = %.4e" % (oA, oB))

    def block(Pi, Pj):
        C = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                acc += (-np.trace(Pi[t] @ tau[t, s] @ Pj[s] @ tau[s, t])).real
            C[dt] = acc / ns
        return C
    CAA = block(PA, PA)
    CBB = block(PB, PB)
    CAB = block(PA, PB)

    def eff(C):
        with np.errstate(all="ignore"):
            return np.log(np.abs(C[:-1]) / np.abs(C[1:]))
    # (b) collinearity: normalized cross-correlation at t0
    cos0 = CAB[T0] / np.sqrt(abs(CAA[T0] * CBB[T0]))
    print("# independence: cos(O_A,O_B) at t0=%d = %.4f  (|1| = collinear)" % (T0, cos0))

    # (c) {O_A, O_B} connected 2x2 GEVP (subtract vacuum constants)
    oo = np.array([[oA * oA, oA * oB], [oA * oB, oB * oB]])
    c = np.zeros((twin, 2, 2))
    for t in range(twin):
        c[t] = np.array([[CAA[t], CAB[t]], [CAB[t], CBB[t]]]) - oo
    C0 = 0.5 * (c[T0] + c[T0].T)
    wv, Uv = np.linalg.eigh(C0)
    keep = wv > 1e-11 * abs(wv).max()             # prune non-positive C0 directions (robust)
    Uk = Uv[:, keep] / np.sqrt(wv[keep])
    ng = int(keep.sum())
    print("# C0 eigenvalues = %s  (kept %d positive)" % (np.array2string(wv, precision=3), ng))
    lam = np.full((twin, ng), np.nan)
    for t in range(twin):
        M = Uk.T @ (0.5 * (c[t] + c[t].T)) @ Uk
        try:
            lam[t] = np.sort(np.linalg.eigvals(M).real)[::-1]
        except Exception:
            pass
    with np.errstate(all="ignore"):
        emg = np.log(lam[:-1] / lam[1:])
    emg = np.pad(emg, ((0, 0), (0, 2 - ng)), constant_values=np.nan)

    eA, eB = eff(CAA), eff(CBB)
    print("\n  dt |  m(O_A) 2pt   m(O_B) 2pt | GEVP lvl0   lvl1")
    for dt in range(1, twin - 2):
        print("  %2d |   %6.3f       %6.3f   |  %6.3f    %6.3f" % (dt, eA[dt], eB[dt], emg[dt, 0], emg[dt, 1]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dts = np.arange(1, twin)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for C, col, mk, lab in [(CAA, "tab:green", "o", r"$\langle O_A O_A\rangle$"),
                            (CBB, "tab:orange", "^", r"$\langle O_B O_B\rangle$")]:
        e = eff(C)
        g = np.isfinite(e) & (C[:-1] * C[1:] > 0)
        ax.plot(dts[:len(e)][g], e[g], color=col, marker=mk, ms=3, lw=1, label=lab)
    g = np.isfinite(emg[:, 0]) & (lam[:-1, 0] * lam[1:, 0] > 0)
    ax.plot(dts[:len(emg)][g], emg[g, 0], color="tab:blue", marker="s", ms=4, lw=1.5,
            label=r"{$O_A,O_B$} GEVP level 0")
    ax.axhline(msig, color="gray", ls="--", lw=1, alpha=0.6, label=r"$2E_0=%.3f$" % msig)
    ax.set_ylim(0.2, 1.1)
    ax.set_xlabel("dt"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d  $(2,2)$ from $O_A$, $O_B$, and their GEVP" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/o_b_test_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
