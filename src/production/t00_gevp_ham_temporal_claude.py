#!/usr/bin/env python3
# t00_gevp_ham_temporal_claude.py
# 2-operator block-Hankel GEVP combining the two fermionic T_00 interpolators (free-limit, 1 cfg):
#   O_H : Hamiltonian naive e.sigma vertex W (r=0),           stencil {(0,1)}       (t00_wilson_kernel build_W)
#   O_T : temporal   M=diag(w)sigma_3, D_t derivative,        stencil {(+1,1/2),(-1,-1/2)}
# Full 2x2 connected cross-correlator (both terms; G=AblkS, Gt=-AblkS off-diag / 1/2 I - AblkS eq-time):
#   C_ij(t,t0) = -sum_{sh_i,sh_j} c_i c_j { Tr[V_i G(t+sh_i,t0) V_j G(t0+sh_j,t)]
#                                          + Tr[V_i^H Gt(t+sh_i,t0) V_j^H Gt(t0+sh_j,t)] }.
# Then the FROZEN block-Hankel core (effmass_axial_tp_l3_perm_hankel) with off0-3 reb1@3 T0=2 (NKEEP=1).
# Run:  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_gevp_ham_temporal_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "final/analysis_axial")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as fg
import t00_wilson_kernel_claude as wk
import t00_stress_temporal_claude as tt
import effmass_axial_tp_l3_perm_hankel_claude as hk

NS = dc.NS
DTMAX = int(os.environ.get("DTMAX", "24"))
OFFS = [int(x) for x in os.environ.get("OFFSETS", "0,3").split(",")]
REBT = int(os.environ.get("REBT", "3"))
NKEEP = int(os.environ.get("NKEEP", "1"))
T0 = int(os.environ.get("T0", "2"))
REF = int(os.environ.get("REF", "5"))               # reference dt for sign/scale normalization
ROVER = 1.0 / 0.189


def corr_matrix(k, ops):
    # ops: list of (V (N,N), stencil [(sh,c),...]).  Returns C[nop,nop,DTMAX] (single loop, both terms).
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(k)
    N = NS * nsite
    Iden = np.eye(N)

    def G(a, b):
        return AblkS(a, b).reshape(N, N)

    def Gt(a, b):
        g = -AblkS(a, b).reshape(N, N)
        if a == b:
            g = 0.5 * Iden - AblkS(a, b).reshape(N, N)
        return g

    nop = len(ops)
    C = np.full((nop, nop, DTMAX), np.nan)
    smax = max(abs(sh) for V, st in ops for sh, c in st)
    for dt in range(DTMAX):
        s_lo = smax
        s_hi = twin - 1 - dt - smax
        if s_hi < s_lo:
            continue
        acc = np.zeros((nop, nop))
        cnt = 0
        for s in range(s_lo, s_hi + 1):
            t0 = s
            t = s + dt
            for i in range(nop):
                Vi, sti = ops[i]
                ViH = Vi.conj().T
                for j in range(nop):
                    Vj, stj = ops[j]
                    VjH = Vj.conj().T
                    a = 0.0
                    b = 0.0
                    for sh_i, c_i in sti:
                        for sh_j, c_j in stj:
                            a += c_i * c_j * np.trace(Vi @ G(t + sh_i, t0) @ Vj @ G(t0 + sh_j, t))
                            b += c_i * c_j * np.trace(ViH @ Gt(t + sh_i, t0) @ VjH @ Gt(t0 + sh_j, t))
                    acc[i, j] += (-(a + b)).real
            cnt += 1
        C[:, :, dt] = acc / cnt
    return C


def main():
    tag = dc.ENS
    dual = dc.dual_areas_from_mesh()
    W, _, _ = wk.build_W("../../geometry/data/", 1, r=0.0)
    M = tt.build_M(dual)
    ops = [(W, [(0, 1.0)]), (M, [(1, 0.5), (-1, -0.5)])]
    ks = dc.KS
    print("# ENS=%s ops={O_H(e.sigma), O_T(temporal)}  off%s reb%d@%d T0=%d  ncfg=%d"
          % (tag, "-".join(map(str, OFFS)), NKEEP, REBT, T0, len(ks)))

    C = np.mean([corr_matrix(k, ops) for k in ks], 0)     # (2,2,DTMAX)
    C = 0.5 * (C + np.swapaxes(C, 0, 1))                   # symmetrize
    nop = C.shape[0]

    # sign+scale normalization: D C D with D=diag(s_i), s_i = sign(C_ii[REF])/sqrt(|C_ii[REF]|) -> C_ii[REF]=1
    s = np.array([np.sign(C[i, i, REF]) / np.sqrt(abs(C[i, i, REF])) for i in range(nop)])
    Cn = C * np.outer(s, s)[:, :, None]
    print("# normalized C_ij[REF=%d]:\n%s" % (REF, np.array2string(Cn[:, :, REF], precision=4)))
    print("# off-diag overlap C01/sqrt(C00 C11) at a few dt:")
    for dt in (3, 5, 8, 11):
        r = Cn[0, 1, dt] / np.sqrt(Cn[0, 0, dt] * Cn[1, 1, dt]) if Cn[0, 0, dt] * Cn[1, 1, dt] > 0 else np.nan
        print("#   dt=%2d  rho=%.4f" % (dt, r))

    # frozen block-Hankel GEVP
    Cts = np.moveaxis(Cn, 2, 0)                            # (DTMAX, nop, nop)
    Big = hk.hankel_off(Cts, OFFS)
    Vop = hk.rebase_vectors(Big, REBT, T0, NKEEP)
    em = hk.rebased_effmass_fixed(Big, Vop, T0, NKEEP)    # (tmax-1, NKEEP)

    # single-operator Hankel effmasses for overlay (each op alone)
    emH = hk.hankel_effmass_scalar(Cn[0, 0], OFFS, REBT, NKEEP, T0, 0.2)[0]
    emT = hk.hankel_effmass_scalar(Cn[1, 1], OFFS, REBT, NKEEP, T0, 0.2)[0]

    print("\n#  t   GEVP(2op)   H-only    T-only   [T00=%.3f]" % (3.0 / ROVER))
    for t in range(em.shape[0]):
        if np.isfinite(em[t, 0]):
            print("  %2d  %8.4f  %8.4f %8.4f" % (t, em[t, 0],
                  emH[t, 0] if t < emH.shape[0] else np.nan,
                  emT[t, 0] if t < emT.shape[0] else np.nan))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tt_ = np.arange(em.shape[0])
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    for val, lab, col in ((2.0 / ROVER, r"$\sigma=2/R$", "tab:green"),
                          (3.0 / ROVER, r"$T_{00}=3/R=0.567$", "tab:red"),
                          (4.0 / ROVER, r"$2m=4/R$", "gray")):
        ax.axhline(val, color=col, ls="--", lw=1, alpha=0.7)
        ax.text(0.2, val + 0.006, lab, fontsize=9, color=col)
    ax.plot(tt_[np.isfinite(emH[:, 0])], emH[np.isfinite(emH[:, 0]), 0], "^:", color="tab:orange", ms=4,
            lw=0.9, label="H-only (Hankel)")
    ax.plot(tt_[np.isfinite(emT[:, 0])], emT[np.isfinite(emT[:, 0]), 0], "v:", color="tab:blue", ms=4,
            lw=0.9, label="T-only (Hankel)")
    ax.plot(tt_[np.isfinite(em[:, 0])], em[np.isfinite(em[:, 0]), 0], "o-", color="tab:red", ms=5.5, lw=1.3,
            label="2-op GEVP")
    ax.set_ylim(0.3, 1.0)
    ax.set_xlim(0, em.shape[0])
    ax.set_xlabel(r"$dt$")
    ax.set_ylabel(r"$a_t m_\mathrm{eff}$")
    ax.set_title(r"$\{O_H,\,O_T\}$ 2-op block-Hankel GEVP  off%s reb%d@%d T0=%d  %s 1cfg"
                 % ("-".join(map(str, OFFS)), NKEEP, REBT, T0, tag), fontsize=10.5)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs("figs", exist_ok=True)
    out = "figs/t00_gevp_ham_temporal_%s_claude.png" % tag
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("\n# -> %s" % out)


if __name__ == "__main__":
    main()
