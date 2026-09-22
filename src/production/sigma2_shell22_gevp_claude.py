#!/usr/bin/env python3
# sigma2_shell22_gevp_claude.py -- single-meson SHELL GEVP to resolve the (2,2) at interacting L2, where the
#   bare ell=3/2 shell operator leaks to m_PS (degeneracy lifts, ranks 5-12 overlap the ground shell).  Build one
#   operator per |eig(M)| rank-window (shell) and GEVP them: the GEVP assigns 2E_0=m_PS to the ell=1/2 op and
#   2E_{3/2}=(2,2) to the ell=3/2 op.  M(t)=tau(t,t)-1/2 I_Nv (stored, anti-hermitian).  Kernel {1,1,1,1}.
#     Op X = P_X(t) = projector onto |eig(M)| ranks in window X.  C_XY(dt)=avg_s -Tr[P_X[t] tau[t,s] P_Y[s] tau[s,t]].
#   Default windows (Nv=24): ell1/2=[0,4], ell3/2=[4,12], ell5/2=[12,24].
#   Run: ENS=<L2> LREF=2 NVDIR=distill_Nv24_v2 WINDOWS=0-4,4-12,12-24 REBT=4 NKEEP=3 T0=2 BINSIZE=10 python3 sigma2_shell22_gevp_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import hankel_rebase_scan_claude as hs

WINDOWS = [tuple(int(x) for x in w.split("-")) for w in os.environ.get("WINDOWS", "0-4,4-12,12-24").split(",")]
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
REBT = int(os.environ.get("REBT", "4"))
NKEEP = int(os.environ.get("NKEEP", str(len(WINDOWS))))
T0 = int(os.environ.get("T0", "2"))
OFFSETS = [int(x) for x in os.environ.get("OFFSETS", "0").split(",")]
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))


def shell_projectors(M):
    w, U = np.linalg.eigh(1j * M)
    order = np.argsort(-np.abs(w))                     # |eig| descending
    Ps = []
    for (lo, hi) in WINDOWS:
        Us = U[:, order[lo:hi]]
        Ps.append(Us @ Us.conj().T)
    return Ps


def one_config(k):
    V, windows = dc.load_peram_windows(k)
    nop = len(WINDOWS)
    C = np.zeros((nop, nop, DTMAX))
    nwin = 0
    for (tsrc0, tau, taugw) in windows:
        twin = tau.shape[0]
        Iv = np.eye(tau.shape[-1])
        P = [shell_projectors(tau[a, a] - 0.5 * Iv) for a in range(twin)]   # P[a][op]
        for dt in range(DTMAX):
            s0s = [s for s in range(twin) if s + dt < twin]
            if not s0s:
                continue
            for a in range(nop):
                for b in range(nop):
                    acc = 0.0
                    for s in s0s:
                        t = s + dt
                        acc += -np.trace(P[t][a] @ tau[t, s] @ P[s][b] @ tau[s, t]).real
                    C[a, b, dt] += acc / len(s0s)
        nwin += 1
    return C / nwin


def gevp(Cmat, Vfix):
    Cts = np.transpose(Cmat, (2, 0, 1))
    Big = hs.hankel_off(Cts, OFFSETS)
    V = hs.staged_project(Big, [(REBT, NKEEP)], T0) if Vfix is None else Vfix
    return hs.rebased_effmass_fixed(Big, V, T0), V


def main():
    tag = dc.ENS.split("nu0")[0]
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# ENS=%s L=%d  SHELL GEVP windows=%s  %d cfg  reb%d@%d T0=%d  m_PS=%.4f 2m_PS=%.4f"
          % (tag, dc.L, WINDOWS, len(ks), NKEEP, REBT, T0, mps, m2ps))
    allC = np.array([one_config(k) for k in ks])
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    em_c, Vfix = gevp(blk.mean(0), None)
    ems = np.array([gevp(np.delete(blk, i, 0).mean(0), Vfix)[0] for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))
    tmax = em_c.shape[0]

    print("\n#  t |  " + "  ".join("m%d(err)   " % n for n in range(NKEEP)))
    for t in range(T0, min(tmax, 16)):
        print("#  %2d | %s" % (t, "  ".join("%7.4f(%.4f)" % (em_c[t, n], em_e[t, n]) for n in range(NKEEP))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(tmax)
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.8)
    ax.text(TMAXPLOT * 0.7, m2ps + 0.012, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="firebrick", ls=":", lw=1.2, alpha=0.8)
    ax.text(TMAXPLOT * 0.7, mps + 0.012, r"$m_{PS}=%.4f$ ($\ell$1/2)" % mps, fontsize=9, color="firebrick")
    cols = ["firebrick", "tab:purple", "tab:blue", "tab:green"]
    mkr = ["o", "D", "s", "^"]
    lab = [r"state 0 ($\ell$1/2 $=m_{PS}$)", r"state 1 ($\ell$3/2 $=(2,2)$)", "state 2 ($\ell$5/2)", "state 3"]
    for n in range(NKEEP):
        g = np.isfinite(em_c[:, n]) & np.isfinite(em_e[:, n]) & (em_e[:, n] < 0.3)
        ax.errorbar(ts[g], em_c[g, n], yerr=em_e[g, n], color=cols[n % 4], marker=mkr[n % 4], ms=6, lw=1.2,
                    capsize=3, label=lab[n] if n < len(lab) else "state %d" % n)
    ax.set_ylim(0.2, 1.1)
    ax.set_xlim(T0, TMAXPLOT)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"single-meson SHELL GEVP (M-eigenspaces)  %s L%d %dcfg  windows=%s"
                 % (tag, dc.L, ncfg, WINDOWS), fontsize=9)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_shell22_gevp_%s_L%d_claude.png" % (tag, dc.L)
    fig.savefig(out, dpi=140)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
