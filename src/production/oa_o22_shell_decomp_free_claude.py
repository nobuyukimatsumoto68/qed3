#!/usr/bin/env python3
# oa_o22_shell_decomp_free_claude.py  [(a) continuum-free mechanism]
# Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 oa_o22_shell_decomp_free_claude.py
#
# Decompose the FREE bilinear loops C_AA (O_A) and C_QQ (O_22) into per-shell-pair (lambda,lambda')
# contributions, by sandwiching the two propagators of the loop with the shell projectors Pi_lambda:
#     C^{K}_{ll'}(dt) = mean_s  -Tr[ P_K(t) (Pi_l tau(t,s) Pi_l) P_K(s) (Pi_l' tau(s,t) Pi_l') ],
#     sum_{l,l'} C^{K}_{ll'} = C_KK.
# A pair (l,l') decays as exp[-(E_l+E_l') t].  This exposes the mechanism:
#   O_22 (K=Q2 tilde_tau)  -> essentially only (2,2) -> clean e^{-2E_2 t} (the genuine (1,1,1,1)=2E_1).
#   O_A  (K=tilde_tau)     -> a mix of shells; the reading is that lambda=1 either decouples or the
#                            mixed-sign shell pairs cancel at late t -> effmass decays to null.
# Free lattice L=1 IS the free theory (no gauge); this is the continuum mechanism via the exact modes.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

DT_DECOMP = int(os.environ.get("DT_DECOMP", "3"))


def all_shell_projectors(tau, twin, Nv, dt_dec):
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
    Pis = []
    Esh = []
    for (i0, j0) in clusters:
        sel = list(range(i0, j0 + 1))
        Pis.append(R[:, sel] @ Rinv[sel, :])
        Esh.append(float(np.mean(E[i0:j0 + 1])))
    return Pis, clusters, Esh


def effmass(C):
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s FREE L=%d  shell-pair decomposition of C_AA (O_A) and C_QQ (O_22)" % (tag, dc.L))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]

    Pis, clusters, Esh = all_shell_projectors(tau, twin, Nv, DT_DECOMP)
    nsh = len(Pis)
    degs = [c[1] - c[0] + 1 for c in clusters]
    print("# shells: degeneracies %s (expect 4,8,12), E_shell %s" % (degs, ["%.4f" % e for e in Esh]))
    Q2 = Pis[1]                                                     # lambda=2 projector
    PQ = [Q2 @ tt[a, a] for a in range(twin)]

    # pre-project the perambulator on each shell for speed: tPi[l][t][s] not stored; do inline
    # per-(l,l') decomposition of both C_AA and C_QQ
    use = min(nsh, 3)
    CAA_ll = np.zeros((use, use, twin))
    CQQ_ll = np.zeros((use, use, twin))
    for dt in range(twin):
        ns = twin - dt
        accA = np.zeros((use, use))
        accQ = np.zeros((use, use))
        for s in range(ns):
            t = s + dt
            # shell-projected forward/backward propagators
            fwd = [Pis[l] @ tau[t, s] @ Pis[l] for l in range(use)]
            bwd = [Pis[l] @ tau[s, t] @ Pis[l] for l in range(use)]
            for l in range(use):
                PAt_f = PA[t] @ fwd[l]
                PQt_f = PQ[t] @ fwd[l]
                for lp in range(use):
                    accA[l, lp] += (-np.trace(PAt_f @ PA[s] @ bwd[lp])).real
                    accQ[l, lp] += (-np.trace(PQt_f @ PQ[s] @ bwd[lp])).real
        CAA_ll[:, :, dt] = accA / ns
        CQQ_ll[:, :, dt] = accQ / ns

    CAA = CAA_ll.sum(axis=(0, 1))
    CQQ = CQQ_ll.sum(axis=(0, 1))

    # report the pair breakdown at a mid time, and the expected pair energies E_l+E_l'
    tm = min(twin - 2, 10)
    print("\n# pair energies E_l + E_l' (l,l' in 1..%d):" % use)
    for l in range(use):
        row = "  ".join("%.3f" % (Esh[l] + Esh[lp]) for lp in range(use))
        print("#   l=%d: %s" % (l + 1, row))
    print("\n# C_AA(t=%d) pair breakdown (fraction of total %.3e):" % (tm, CAA[tm]))
    for l in range(use):
        row = "  ".join("%+8.3f" % (CAA_ll[l, lp, tm] / (CAA[tm] + 1e-300)) for lp in range(use))
        print("#   l=%d: %s" % (l + 1, row))
    print("# C_QQ(t=%d) pair breakdown (fraction of total %.3e):" % (tm, CQQ[tm]))
    for l in range(use):
        row = "  ".join("%+8.3f" % (CQQ_ll[l, lp, tm] / (CQQ[tm] + 1e-300)) for lp in range(use))
        print("#   l=%d: %s" % (l + 1, row))

    # effmass of total and of the dominant single pairs
    emA = effmass(CAA)
    emQ = effmass(CQQ)
    emA_11 = effmass(CAA_ll[0, 0])
    emA_22 = effmass(CAA_ll[1, 1])
    emQ_22 = effmass(CQQ_ll[1, 1])
    print("\n#  t |  m_A(tot)  m_A(1,1)  m_A(2,2) |  m_Q(tot)  m_Q(2,2) | 2E_1=%.3f 2E_2=%.3f" %
          (2 * Esh[0], 2 * Esh[1]))
    for t in range(1, min(twin - 1, 26)):
        def fmt(em, C):
            return "%8.4f" % em[t] if np.isfinite(em[t]) and C[t] * C[t + 1] > 0 else "  ----  "
        print("# %2d | %s  %s  %s | %s  %s" %
              (t, fmt(emA, CAA), fmt(emA_11, CAA_ll[0, 0]), fmt(emA_22, CAA_ll[1, 1]),
               fmt(emQ, CQQ), fmt(emQ_22, CQQ_ll[1, 1])))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(twin)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    # left: |contributions| of shell pairs to C_AA
    for l in range(use):
        for lp in range(l, use):
            lab = "(%d,%d)" % (l + 1, lp + 1)
            c = CAA_ll[l, lp] + (CAA_ll[lp, l] if lp != l else 0.0)
            g = np.abs(c) > 0
            ax[0].semilogy(ts[g], np.abs(c[g]), marker=".", ms=3, lw=1, label=lab)
    ax[0].semilogy(ts, np.abs(CAA), "k-", lw=1.5, label="total")
    ax[0].set_title(r"$C_{AA}$ ($O_A$) shell-pair $|$contrib$|$")
    ax[0].set_xlabel("t")
    ax[0].legend(fontsize=8, ncol=2)
    # right: effective masses
    gA = np.isfinite(emA) & (CAA[:-1] * CAA[1:] > 0)
    gQ = np.isfinite(emQ) & (CQQ[:-1] * CQQ[1:] > 0)
    ax[1].plot(np.arange(twin - 1)[gA], emA[gA], "r-o", ms=3, lw=1, label=r"$O_A$ total")
    ax[1].plot(np.arange(twin - 1)[gQ], emQ[gQ], "b-s", ms=3, lw=1, label=r"$O_{22}$ total")
    # physical meson-scale references (perambulator K(dt) energies are in a different internal unit)
    for y, lab in [(0.378, r"$2m_\sigma^{(\lambda1)}$"), (0.52, r"$(2,2){=}2E_1$"), (0.756, r"$2m_\sigma$")]:
        ax[1].axhline(y, color="k", ls=":", lw=0.8, alpha=0.4)
        ax[1].text(twin - 3, y + 0.008, lab, fontsize=8, alpha=0.7)
    ax[1].set_ylim(-0.05, 1.0)
    ax[1].set_xlabel("t")
    ax[1].set_ylabel(r"$m_\mathrm{eff}$")
    ax[1].set_title("effective masses")
    ax[1].legend(fontsize=9)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/oa_o22_shell_decomp_free_L%d_claude.png" % dc.L
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
