#!/usr/bin/env python3
# sigma2_shell22_claude.py -- the (2,2) single-meson DIAGONAL correlator/effmass via the ell=3/2 SHELL PROJECTOR
#   of the stored anti-hermitian M(t) = tau(t,t) - 1/2 I_Nv (mode space).  Kernel + selection from the {1,1,1,1}
#   agent (verified free L1 -> 2E_{3/2}=0.556, no ground leak).  NO new solve, NO Xi -- built from /peram/tau.
#     Phi22(t) = P_shell(t) = projector onto |eig(M)|-ranks [LO,HI) (default 4..12 = the 8-mode ell=3/2 shell).
#     C(dt) = translation-avg  -Tr[ Phi22[t] tau[t,s] Phi22[s] tau[s,t] ] ,  t=s+dt   (all Nv x Nv, mode space).
#   Run: ENS=<L2> LREF=2 NVDIR=distill_Nv24_v2 LO=4 HI=12 BINSIZE=10 python3 sigma2_shell22_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

LO = int(os.environ.get("LO", "4"))          # rank window [LO,HI) on |eig(M)| descending; 4..12 = ell3/2 (2,2)
HI = int(os.environ.get("HI", "12"))
DTMAX = int(os.environ.get("DTMAX", "20"))
BINSIZE = int(os.environ.get("BINSIZE", "10"))
NCFG = int(os.environ.get("NCFG", "0"))
TMAXPLOT = int(os.environ.get("TMAXPLOT", "18"))


def shell_projector(M):
    w, U = np.linalg.eigh(1j * M)                     # iM hermitian; w real, U orthonormal
    sel = np.argsort(-np.abs(w))[LO:HI]               # |eig| ranks LO..HI = the shell
    Us = U[:, sel]
    return Us @ Us.conj().T                            # Nv x Nv orthogonal projector


def one_config(k):
    V, windows = dc.load_peram_windows(k)
    C = np.zeros(DTMAX)
    nwin = 0
    for (tsrc0, tau, taugw) in windows:
        twin = tau.shape[0]
        Nv = tau.shape[-1]
        Iv = np.eye(Nv)
        Phi = [shell_projector(tau[a, a] - 0.5 * Iv) for a in range(twin)]
        for dt in range(DTMAX):
            s0s = [s for s in range(twin) if s + dt < twin]
            if not s0s:
                continue
            acc = 0.0
            for s in s0s:
                t = s + dt
                acc += -np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t]).real
            C[dt] += acc / len(s0s)
        nwin += 1
    return C / nwin


def main():
    tag = dc.ENS.split("nu0")[0]
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    m2ps = float(os.environ.get("M2PS", "0.7054" if dc.L == 2 else "0.6418"))
    print("# ENS=%s L=%d NVDIR=%s  (2,2) shell projector ranks [%d,%d)  %d cfg  m_PS=%.4f 2m_PS=%.4f"
          % (tag, dc.L, os.environ.get("NVDIR"), LO, HI, len(ks), mps, m2ps))
    allC = np.array([one_config(k) for k in ks])
    ncfg = allC.shape[0]
    nb = ncfg // BINSIZE
    blk = np.array([allC[i * BINSIZE:(i + 1) * BINSIZE].mean(0) for i in range(nb)])
    C = blk.mean(0)
    with np.errstate(all="ignore"):
        em_c = np.log(C[:-1] / C[1:])
        ems = np.array([np.log(np.delete(blk, i, 0).mean(0)[:-1] / np.delete(blk, i, 0).mean(0)[1:]) for i in range(nb)])
    em_e = np.sqrt((nb - 1) * np.mean((ems - ems.mean(0)) ** 2, axis=0))

    print("\n#  t |  a_t m_eff(err)   [ (2,2) = 2E_{3/2}; free 0.556 ]")
    for t in range(2, min(len(em_c), 16)):
        print("#  %2d | %7.4f(%.4f)" % (t, em_c[t], em_e[t]))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ts = np.arange(len(em_c))
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    ax.axhline(m2ps, color="gray", ls="--", lw=1.1, alpha=0.8)
    ax.text(TMAXPLOT * 0.7, m2ps + 0.012, r"$2m_{PS}=%.4f$" % m2ps, fontsize=9, color="gray")
    ax.axhline(mps, color="dimgray", ls=":", lw=1.1, alpha=0.8)
    ax.text(TMAXPLOT * 0.7, mps + 0.012, r"$m_{PS}=%.4f$" % mps, fontsize=9, color="dimgray")
    g = np.isfinite(em_c) & np.isfinite(em_e) & (em_e < 0.3)
    ax.errorbar(ts[g], em_c[g], yerr=em_e[g], color="tab:purple", marker="D", ms=6, lw=1.2, capsize=3,
                label=r"$(2,2)$ shell ($\ell{=}3/2$, ranks %d-%d)" % (LO + 1, HI))
    ax.set_ylim(0.2, 1.1)
    ax.set_xlim(1, TMAXPLOT)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"(2,2) single-meson shell-projector effmass  %s L%d  %dcfg  ranks[%d,%d)"
                 % (tag, dc.L, ncfg, LO, HI), fontsize=10)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/sigma2_shell22_%s_L%d_ranks%d-%d_claude.png" % (tag, dc.L, LO, HI)
    fig.savefig(out, dpi=140)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
