#!/usr/bin/env python3
# sigma2_shell_Uscan_claude.py -- scan the SINGLE-EIGENVECTOR kernels K_i = u_i u_i^H one by one (u_i = distillation
#   eigenvectors, NO M anywhere) and show the overlap patterns.  In the U basis the full Nv x Nv correlator matrix is
#     C_ij(t,s) = -Re[ tau(t,s)_ij tau(s,t)_ji ]      (elementwise product).
#   Outputs: (1) normalized matrix R_ij = C_ij/sqrt(C_ii C_jj) at dt=DTN (block pattern = natural slicing);
#            (2) diagonal effmass per eigenvector at a few dt; (3) ground-state overlap fraction per eigenvector,
#                f_i = C_ii(TL) / max_j C_jj(TL)  (kernels are unit-norm; dt=0 is pure contact, not used).
#   Run: ENS=<L1> LREF=1 NVDIR=distill_Nv24_sym NCFG=200 python3 sigma2_shell_Uscan_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as G

DTMAX = int(os.environ.get("DTMAX", "20"))
NCFG = int(os.environ.get("NCFG", "0"))
DTN = int(os.environ.get("DTN", "4"))
TL = int(os.environ.get("TL", "12"))
EFFT = [int(x) for x in os.environ.get("EFFT", "2,5,9").split(",")]


def one_config(k):
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(k, 0)
    nv = tau.shape[-1]
    C = np.zeros((nv, nv, DTMAX))
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        acc = np.zeros((nv, nv))
        for s in s0s:
            t = s + dt
            acc += -(tau[t, s] * tau[s, t].T).real
        C[:, :, dt] = acc / len(s0s)
    return C


def main():
    tag = dc.ENS.split("nu0")[0]
    ks = dc.KS if NCFG == 0 else dc.KS[:NCFG]
    mps = float(os.environ.get("M_PS", "0.3527" if dc.L == 2 else "0.3209"))
    print("# ENS=%s L=%d NVDIR=%s  SINGLE-EVEC kernel scan  %d cfg" % (tag, dc.L, dc.NVDIR, len(ks)))
    os.makedirs("cache_claude", exist_ok=True)
    cf = "cache_claude/Uscan_%s_L%d_%s_n%d.npy" % (tag, dc.L, dc.NVDIR, len(ks))
    if os.path.exists(cf):
        allC = np.load(cf)
    else:
        allC = np.array([one_config(k) for k in ks])
        np.save(cf, allC)
    C = allC.mean(0)
    nv = C.shape[0]
    d = np.array([C[i, i] for i in range(nv)])
    R = C[:, :, DTN] / np.sqrt(np.abs(np.outer(d[:, DTN], d[:, DTN])))
    with np.errstate(all="ignore"):
        em = np.log(d[:, :-1] / d[:, 1:])
    f = d[:, TL] / np.max(d[:, TL])

    print("\n# per-eigenvector: index | C_ii(1) | effmass at t=%s | relative large-t amplitude f_i (TL=%d)" % (EFFT, TL))
    for i in range(nv):
        print("#  %2d | %10.4e | %s | %8.4f" % (i, d[i, 1], "  ".join("%7.4f" % em[i, t] for t in EFFT), f[i]))
    print("\n# normalized matrix R_ij at dt=%d (rows i, cols j; x100, rounded):" % DTN)
    for i in range(nv):
        print("#  %2d | %s" % (i, " ".join("%4d" % int(round(100 * R[i, j])) for j in range(nv))))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs("figs", exist_ok=True)
    stem = "figs/sigma2_shell_Uscan_%s_%s_L%d_%s_claude.png"

    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    im = ax.imshow(R, cmap="RdBu_r", vmin=-1, vmax=1, origin="upper")
    fig.colorbar(im, ax=ax, label=r"$C_{ij}/\sqrt{C_{ii}C_{jj}}$ at $dt=%d$" % DTN)
    ax.set_xlabel("eigenvector index j")
    ax.set_ylabel("eigenvector index i")
    ax.set_title("single-evec kernel correlation pattern  %s L%d %s" % (tag, dc.L, dc.NVDIR), fontsize=9)
    fig.tight_layout()
    out = stem % ("R", tag, dc.L, dc.NVDIR)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("# -> %s" % out)

    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    cols = ["firebrick", "tab:blue", "tab:green", "tab:purple"]
    mkr = ["o", "s", "^", "D"]
    for n, t in enumerate(EFFT):
        ax.plot(np.arange(nv), em[:, t], color=cols[n % 4], marker=mkr[n % 4], ms=5, lw=1.0, label="t=%d" % t)
    ax.axhline(mps, color="k", ls=":", lw=1.0, alpha=0.6)
    ax.set_xlabel("eigenvector index i")
    ax.set_ylabel(r"diagonal $a_t m_\mathrm{eff}$")
    ax.set_title("single-evec kernel diagonal effmass  %s L%d %s" % (tag, dc.L, dc.NVDIR), fontsize=9)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = stem % ("effm", tag, dc.L, dc.NVDIR)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("# -> %s" % out)

    fig, ax = plt.subplots(figsize=(9.0, 5.6))
    ax.plot(np.arange(nv), f, color="firebrick", marker="o", ms=5, lw=1.0)
    ax.set_xlabel("eigenvector index i")
    ax.set_ylabel(r"relative large-$t$ amplitude $C_{ii}(%d)/\max_j C_{jj}(%d)$" % (TL, TL))
    ax.set_title("single-evec kernel overlap with $m_{PS}$  %s L%d %s" % (tag, dc.L, dc.NVDIR), fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = stem % ("gsfrac", tag, dc.L, dc.NVDIR)
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("# -> %s" % out)


if __name__ == "__main__":
    main()
