#!/usr/bin/env python3
# sigma2_cross_asymm_gevp_claude.py
#   ASYMMETRIC (non-hermitian) GEVP applied DIRECTLY to the PP-FF cross block of the sigma^2 correlator.
#
#   Cross-correlator (O_i = PP-geom, O~_j = FF-geom, both hermitian, but O_i != O~_j so C is NOT symmetric):
#       C_{i j}(t) = <0| O_i(t) O~_j(0) |0> = sum_n <0|O_i|n> e^{-E_n t} <n|O~_j|0>* = (g Lambda(t) g~^dag)_{ij}
#   The two-point pencil  C(t0)^{-1} C(t)  has eigenvalues lambda_n(t,t0) = e^{-E_n (t-t0)} with distinct
#   LEFT/RIGHT eigenvectors (asymmetric).  "Politerated" (block-Hankel / GPOF) enlargement uses
#   {C(t), C(t+1), C(t+2)} = offsets [0,1] to add resolving power, exactly as in NM's asymm_gevp note.
#
#   Physics test: if the two-meson P^2 and F^2 are DIFFERENT states (<2m,P|2m,F>=0) mixing only through
#   the single-meson (2,2), then the CROSS channel carries ONLY the single-meson pole -> this GEVP should
#   plateau at the single-meson mass (L1 0.556 / L2 ~0.69), with the two-meson absent.
#
#   Method refs: asymmetric/non-hermitian GEVP for cross-correlators (NM asymm_gevp note); generalized
#   pencil-of-function / block-Hankel  Hua-Sarkar GPOF (IEEE 1990), Aubin-Orginos 1010.0202;
#   GEVP Luscher-Wolff.
#   Run: LTAG=L1 OFFSETS=0,1 T0=3 python3 sigma2_cross_asymm_gevp_claude.py
#        LTAG=L2 CACHE=... M_PS=0.393 STATE_A=0.69 TWO_MESON=0.786 python3 sigma2_cross_asymm_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import scipy.linalg as sla
import hankel_rebase_scan_claude as hs

LTAG = os.environ.get("LTAG", "L1")
CACHE = os.environ.get("CACHE",
    "sigma2_flavor_cache_claude/sigma2_flavorgeom_FULL_free_1cfg_d1_claude.npy")
T0 = int(os.environ.get("T0", "3"))
NOHANKEL = int(os.environ.get("NOHANKEL", "0"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", "0"))                # 0 = full-spectrum fixed-t0 GEVP; >0 = biorthogonal rebase keep NKEEP
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0,1").split(",")]
STATE_1111 = float(os.environ.get("STATE_A", "0.556"))   # single-meson (2,2): L1 0.556, L2 ~0.69
TWO_MESON = float(os.environ.get("TWO_MESON", "0.756"))  # 2 m_PS: L1 0.756, L2 0.786
NSHOW = int(os.environ.get("NSHOW", "3"))
DTMAX = int(os.environ.get("DTMAX", "16"))


def rebase_biortho(Big, tr, t0, nkeep):
    # asymmetric rebase: leading nkeep RIGHT/LEFT eigenvectors of the pencil (C(tr), C(t0)),
    # biorthonormalized in the C(t0) metric (Vl^H C(t0) Vr = I).  Non-hermitian analog of
    # staged_project's rebase_vectors.
    A = Big[tr]
    B = Big[t0]
    w, vl, vr = sla.eig(A, B, left=True, right=True)
    order = np.argsort(-w.real)                                    # largest lambda = lowest energy first
    idx = [j for j in order if abs(w[j].imag) < 0.2 * abs(w[j].real) + 1e-12 and w[j].real > 0][:nkeep]
    Vr = vr[:, idx]
    Vl = vl[:, idx]
    d = np.diag(Vl.conj().T @ B @ Vr)
    Vl = Vl / np.conj(d)[None, :]                                  # so Vl^H B Vr = I on the diagonal
    return Vl, Vr


def rebased_effmass(Big, Vl, Vr, t0):
    # project onto the rebased subspace, then effmass.  nkeep=1 -> scalar channel log-ratio;
    # nkeep>1 -> asymmetric fixed-t0 GEVP in the reduced space.
    tmax = Big.shape[0]
    nk = Vr.shape[1]
    Cr = np.array([Vl.conj().T @ Big[t] @ Vr for t in range(tmax)])   # (tmax, nk, nk), complex
    if nk == 1:
        c = Cr[:, 0, 0].real
        with np.errstate(all="ignore"):
            return np.log(c[:-1] / c[1:])[:, None]
    E = np.full((tmax, nk), np.nan)
    A0 = Cr[t0]
    for t in range(tmax):
        try:
            lam = sla.eig(np.linalg.solve(A0, Cr[t]), right=False)
            lam = lam[np.abs(lam.imag) < 0.15 * np.abs(lam.real) + 1e-12].real
            lam = np.sort(lam[lam > 1e-8])[::-1]
            E[t, :len(lam)] = lam
        except Exception:
            pass
    with np.errstate(all="ignore"):
        return np.log(E[:-1] / E[1:])


def gevp_energies(Big, t0):
    # fixed-t0 asymmetric GEVP: eig of C(t0)^{-1} C(t); energies E_n(t) = -ln(lambda_n)/(t-t0)
    tmax = Big.shape[0]
    N = Big.shape[1]
    A0 = Big[t0]
    E = np.full((tmax, N), np.nan)
    for t in range(tmax):
        if t == t0 or np.any(~np.isfinite(Big[t])) or np.any(~np.isfinite(A0)):
            continue
        try:
            lam = sla.eig(np.linalg.solve(A0, Big[t]), left=False, right=False)
        except Exception:
            continue
        lam = lam[np.abs(lam.imag) < 0.15 * np.abs(lam.real) + 1e-12]   # keep ~real eigenvalues
        lam = lam.real
        lam = lam[lam > 1e-8]                                           # positive (physical) branch
        en = -np.log(lam) / (t - t0)
        en = np.sort(en)                                               # lowest energy first
        E[t, :len(en)] = en[:N]
    return E


def main():
    C = np.load(CACHE)
    C = 0.5 * (C + np.swapaxes(C, 1, 2))
    M = C.mean(0)                                          # (9,9,twin)  full symmetrized
    # asymmetric cross series: rows = PP-geom (0,1,2), cols = FF-geom (3,4,5)
    Ccross = np.transpose(M[0:3, 3:6, :], (2, 0, 1))       # (twin,3,3), NOT symmetric

    if NOHANKEL:
        Big = Ccross
        htag_m = "no Hankel"
    else:
        Big = hs.hankel_off(Ccross, OFFSETS)               # politerated block-Hankel
        htag_m = "politerated Dt=%s" % OFFSETS

    if NKEEP > 0:
        Vl, Vr = rebase_biortho(Big, REBT, T0, NKEEP)
        E = rebased_effmass(Big, Vl, Vr, T0)
        method = "asymm GEVP %s reb%d@%d T0=%d" % (htag_m, NKEEP, REBT, T0)
    else:
        E = gevp_energies(Big, T0)
        method = "asymm GEVP %s T0=%d" % (htag_m, T0)
    tmax = E.shape[0]
    nshow = min(NSHOW, E.shape[1])

    print("# FREE %s  PP-FF CROSS asymmetric GEVP  %s" % (LTAG, method))
    print("# refs: single-meson (2,2) = %.3f ;  two-meson 2m_PS = %.3f" % (STATE_1111, TWO_MESON))
    print("#  t | " + "  ".join("E%d" % n for n in range(nshow)))
    for t in range(T0 + 1, min(tmax, DTMAX)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f" % E[t, n] if np.isfinite(E[t, n]) else "  ---  "
                                            for n in range(nshow))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.axhline(STATE_1111, color="tab:orange", ls="--", lw=1, alpha=0.7)
    ax.text(DTMAX * 0.45, STATE_1111 + 0.008, r"single-meson (2,2) $=%.3f$" % STATE_1111,
            color="tab:orange", fontsize=9)
    ax.axhline(TWO_MESON, color="tab:blue", ls="--", lw=1, alpha=0.7)
    ax.text(DTMAX * 0.45, TWO_MESON + 0.008, r"two-meson $2m_{PS}=%.3f$" % TWO_MESON,
            color="tab:blue", fontsize=9)
    cols = ["tab:red", "tab:purple", "tab:brown"]
    mkr = ["o", "s", "^"]
    for n in range(nshow):
        g = np.isfinite(E[:, n])
        ax.plot(ts[g], E[g, n], color=cols[n % 3], marker=mkr[n % 3], ms=5, lw=1.1, label="cross state %d" % n)
    ax.set_ylim(0.3, 1.1)
    ax.set_xlim(T0, min(tmax, DTMAX))
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t E_n$")
    ax.set_title(r"FREE %s  PP-FF cross asymmetric GEVP  (%s)" % (LTAG, method))
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    htag = "noHankel" if NOHANKEL else "off" + "".join(str(o) for o in OFFSETS)
    if NKEEP > 0:
        htag += "_reb%d_%d" % (NKEEP, REBT)
    out = "figs/sigma2_cross_asymm_gevp_%s_%s_T0%d_claude.png" % (LTAG, htag, T0)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
