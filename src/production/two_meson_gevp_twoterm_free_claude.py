#!/usr/bin/env python3
# two_meson_gevp_twoterm_free_claude.py  [FREE -- {1, sigma_PS, sigma_PS^2} GEVP, FULL G4 = 1st + 2nd term]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 two_meson_gevp_twoterm_free_claude.py
#
# qed3int_v3-4.pdf Eq (5.5): every PS correlator = <S..>(tau) + <Stilde..>(DH), DH = D_ov^{-dag}, tau_DH=delta I-tau.
# Blocks (both terms): C11 (single), C12 (triangle), C22 = G4 (two-meson).  Identity carries the vacuum
# (one-points from large-t plateaus).  For PS the 2nd term = conj(1st) -> uniform x2 -> eigenvalues unchanged;
# we print the single-term GEVP alongside to confirm.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

T0 = 3


def solve_gevp(Cts, T0, tol=1e-11):
    twin = Cts.shape[0]
    C0 = 0.5 * (Cts[T0] + Cts[T0].T)
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
    print("# ENS=%s  FREE L=%d  {1, sigma_PS, sigma_PS^2} GEVP  FULL G4 (1st+2nd term)  m_sigma~%.3f"
          % (tag, dc.L, msig))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = tau.shape[-1]
    Iv = np.eye(Nv)
    tau_DH = np.conj(tau).transpose(1, 0, 3, 2)         # D^{-dag}
    legs = [tau, tau_DH]
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    def build(leglist):
        C11 = np.zeros(twin); C12 = np.zeros(twin); C22 = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            a11 = a12 = a22 = 0.0
            for s in range(ns):
                t = s + dt
                for lg in leglist:
                    a11 += (-np.trace(Phi[t] @ lg[t, s] @ Phi[s] @ lg[s, t])).real
                    tss = lg[s, s] - 0.5 * Iv                       # tilde on the sigma^2 source loop
                    a12 += (-2.0 * np.trace(Phi[s] @ tss @ Phi[s] @ lg[s, t] @ Phi[t] @ lg[t, s])).real
                    a22 += (dc.W10 * de.diags_pair(Phi, lg, s, t)).sum().real
            C11[dt], C12[dt], C22[dt] = a11 / ns, a12 / ns, a22 / ns
        return C11, C12, C22

    C11, C12, C22 = build(legs)                          # two-term (full G4)
    c11, c12, c22 = build([tau])                         # single-term (1st only), for comparison

    def assemble(C11, C12, C22):
        o1 = 0.0                                         # <sigma_PS> = 0 (contact-subtracted)
        o2 = np.sqrt(max(C22[twin - 6:twin - 1].mean(), 0.0))   # <sigma_PS^2> from plateau
        Cts = np.zeros((twin, 3, 3))
        for dt in range(twin):
            Cts[dt] = np.array([[1.0, o1, o2], [o1, C11[dt], C12[dt]], [o2, C12[dt], C22[dt]]])
        return Cts, o2

    Cts, o2 = assemble(C11, C12, C22)
    cts, o2s = assemble(c11, c12, c22)
    print("# <sigma_PS^2>: full G4 = %.4e , single-term = %.4e  (ratio %.3f ~ sqrt2)" % (o2, o2s, o2 / o2s))
    lam, nl = solve_gevp(Cts, T0)
    lam2, nl2 = solve_gevp(cts, T0)
    with np.errstate(all="ignore"):
        em = np.log(lam[:-1] / lam[1:])
        em2 = np.log(lam2[:-1] / lam2[1:])

    print("\n  t  | FULL G4 {1,PS,PS^2}: m0     m1     m2   | single-term: m0     m1     m2")
    for t in range(T0 + 1, twin - 2):
        r = "  ".join("%6.3f" % em[t, i] if i < nl and np.isfinite(em[t, i]) else "  --- " for i in range(3))
        r2 = "  ".join("%6.3f" % em2[t, i] if i < nl2 and np.isfinite(em2[t, i]) else "  --- " for i in range(3))
        print("  %2d | %s  | %s" % (t, r, r2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(0, twin - 1)
    cols = ["tab:gray", "tab:red", "tab:blue"]
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for i in range(nl):
        g = np.isfinite(em[:, i]) & (lam[:-1, i] * lam[1:, i] > 0)
        ax.plot(ts[g], em[g, i], color=cols[i % 3], marker="o", ms=3, lw=1, label="level %d" % i)
    for y in [0.0, msig, 2 * msig]:
        ax.axhline(y, color="k", ls="--", lw=0.8, alpha=0.4)
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlabel("t"); ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"FREE L=%d {$1,\sigma_{PS},\sigma_{PS}^2$} GEVP, full $G_4$ (1st+2nd)" % dc.L)
    ax.legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/two_meson_gevp_twoterm_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
