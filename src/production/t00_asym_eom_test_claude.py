#!/usr/bin/env python3
# t00_asym_eom_test_claude.py
# Test whether the ANTISYMMETRIZED temporal T_00 (paper Eq IV.4, derivative on BOTH legs, t-even) overlaps the
# spatial Hamiltonian O_H (e.sigma).  My earlier O_T had D_t on xi only -> t-ODD -> orthogonal to O_H.
# General term-list contraction: an operator = list of terms (kind, a, b, c, V) meaning
#   c * barfield^H(t+a) V field(t+b),  kind 'ex' (eta^H..xi, prop G=AblkS) or 'xe' (xi^H..eta, prop Gt=-AblkS).
# Connected cross of term_i (sink) and term_j (source), same kind:
#   'ex':  -c_i c_j Tr[V_i G(t+b_i, t0+a_j) V_j G(t0+b_j, t+a_i)]
#   'xe':  -c_i c_j Tr[V_i Gt(t+b_i, t0+a_j) V_j Gt(t0+b_j, t+a_i)]   (cross kind = 0).
# Operators (M = diag(w) sigma_3, w=A_x Y00; W = r=0 e.sigma kernel):
#   O_H    = eta^H W xi + xi^H W^H eta                                  (static, t-even = T_00)
#   O_T    = eta^H M (D_t xi) + h.c.                                    (D_t on xi only, t-ODD)
#   O_asym = eta^H M (D_t xi) - (D_t eta^H) M xi + h.c.                 (antisymmetrized, t-even, Eq IV.4)
# EOM check: rho_HA = C_HA/sqrt(C_HH C_AA) -> |rho|=1 if O_asym ~ O_H (same state).
# Run:  ENS=free NVDIR=distill_Nv24 LREF=1 python3 t00_asym_eom_test_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import fs_gevp_point_claude as fg
import t00_wilson_kernel_claude as wk
import t00_stress_temporal_claude as tt

NS = dc.NS
DTMAX = int(os.environ.get("DTMAX", "24"))
ROVER = 1.0 / 0.189


def corr_ops(k, A, B):
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

    smax = max(abs(v) for term in (A + B) for v in (term[1], term[2]))
    C = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        s_lo = smax
        s_hi = twin - 1 - dt - smax
        if s_hi < s_lo:
            continue
        acc = 0.0
        cnt = 0
        for s in range(s_lo, s_hi + 1):
            t0 = s
            t = s + dt
            val = 0.0
            for (ki, ai, bi, ci, Vi) in A:
                for (kj, aj, bj, cj, Vj) in B:
                    if ki != kj:
                        continue
                    if ki == "ex":
                        val += -ci * cj * np.trace(Vi @ G(t + bi, t0 + aj) @ Vj @ G(t0 + bj, t + ai))
                    else:
                        val += -ci * cj * np.trace(Vi @ Gt(t + bi, t0 + aj) @ Vj @ Gt(t0 + bj, t + ai))
            acc += val
            cnt += 1
        C[dt] = (acc / cnt).real
    return C


def main():
    dual = dc.dual_areas_from_mesh()
    W, _, _ = wk.build_W("../../geometry/data/", 1, r=0.0)
    WH = W.conj().T
    M = tt.build_M(dual)
    k = dc.KS[0]

    O_H = [("ex", 0, 0, 1.0, W), ("xe", 0, 0, 1.0, WH)]
    O_T = [("ex", 0, 1, 0.5, M), ("ex", 0, -1, -0.5, M),
           ("xe", 0, 1, 0.5, M), ("xe", 0, -1, -0.5, M)]
    # antisymmetrized: P (D_t on xi) minus Q (D_t on eta^H); both ex and xe
    O_A = [("ex", 0, 1, 0.5, M), ("ex", 0, -1, -0.5, M),      # P: eta^H(t) M xi(t+-1)
           ("ex", 1, 0, -0.5, M), ("ex", -1, 0, 0.5, M),       # Q: -(eta^H(t+-1)) M xi(t)
           ("xe", 0, 1, 0.5, M), ("xe", 0, -1, -0.5, M),
           ("xe", 1, 0, -0.5, M), ("xe", -1, 0, 0.5, M)]

    C_HH = corr_ops(k, O_H, O_H)
    C_TT = corr_ops(k, O_T, O_T)
    C_AA = corr_ops(k, O_A, O_A)
    C_HT = corr_ops(k, O_H, O_T)
    C_HA = corr_ops(k, O_H, O_A)

    print("# free 1cfg.  rho_XY = C_XY/sqrt(|C_XX C_YY|)  (|rho|=1 => same state; EOM check for H vs asym)")
    print("#  dt |   C_HH        C_AA        C_HA       rho_HA   |  C_HT      rho_HT")
    for dt in range(2, DTMAX - 2):
        if not (np.isfinite(C_HH[dt]) and np.isfinite(C_AA[dt])):
            continue
        dHA = np.sqrt(abs(C_HH[dt] * C_AA[dt]))
        rHA = C_HA[dt] / dHA if dHA > 0 else np.nan
        dHT = np.sqrt(abs(C_HH[dt] * C_TT[dt]))
        rHT = C_HT[dt] / dHT if dHT > 0 else np.nan
        print("#  %2d | % .3e % .3e % .3e  %+7.4f  | % .1e  %+7.4f"
              % (dt, C_HH[dt], C_AA[dt], C_HA[dt], rHA, C_HT[dt], rHT))

    # effmass of the antisymmetrized operator (sign-normalized)
    sA = np.sign(C_AA[3])
    with np.errstate(all="ignore"):
        emA = np.log((sA * C_AA)[:-1] / (sA * C_AA)[1:])
    print("\n# antisymmetrized T00 effmass (sign %+d):  [T00=%.3f]" % (sA, 3.0 / ROVER))
    for dt in range(2, DTMAX - 2):
        if np.isfinite(emA[dt]):
            print("#   dt=%2d  m=%7.4f" % (dt, emA[dt]))


if __name__ == "__main__":
    main()
