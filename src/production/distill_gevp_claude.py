#!/usr/bin/env python3
# distill_gevp_claude.py
#
# CHUNK 1 of the coupled asymmetric GEVP (gevp_lanczos_impl_plan_claude.md): the sigma\sigma 2x2 GEVP
# {sigma_PS^2, sigma_FS^2}.  Reuses the distillation four-point (distill_contract_claude.py).
#
# For a correlator <O_a(t) O_b(0)> (O_a = sigma_a sigma_a), the four-point = <S_4> + <S~_4> (Eq 5.5):
#   S_4  : all legs tau  (same for every channel).
#   S~_4 : a leg's type is set by its SINK vertex -- vertex at timeslice t (sink operator a) uses a's leg
#          (tau if PS, -tau' if FS); vertex at 0 (source operator b) uses b's.  So compute_diags is
#          generalized to TWO sink-keyed leg-providers (leg_t for sink-at-t legs, leg_0 for sink-at-0).
# Hence C[PS,FS] != C[FS,PS]: the matrix is ASYMMETRIC (sigma_FS is non-Hermitian), so we solve the
# asymmetric GEVP C(t) v = lambda C(t0) v (distinct left/right eigenvectors).
#
# Vacuum subtraction: 0++ has vacuum quantum numbers, so subtract the large-dt plateau (= <O_a><O_b>).

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

CH = ["PS", "FS"]                 # channel index 0=PS, 1=FS


def compute_diags2(Phi, leg_t, leg_0, twin):
    # 10 diagrams with sink-keyed legs: legs whose SINK is at t use leg_t; sink at 0 use leg_0.
    Phi0 = Phi[0]
    t00 = leg_0[0, 0]            # sink 0
    DS0 = np.trace(Phi0 @ t00)
    DpS0 = np.trace(Phi0 @ t00 @ Phi0 @ t00)
    diags = np.zeros((10, twin), complex)
    for dt in range(twin):
        Phit = Phi[dt]
        ttt = leg_t[dt, dt]      # sink t
        t0t = leg_0[0, dt]       # sink 0
        tt0 = leg_t[dt, 0]       # sink t
        DSt = np.trace(Phit @ ttt)
        DpSt = np.trace(Phit @ ttt @ Phit @ ttt)
        M = Phi0 @ t0t @ Phit @ tt0
        CS = np.trace(M)
        TS = np.trace(M @ M)
        VS_0t = np.trace(Phi0 @ t00 @ M)
        SS_0t = np.trace(Phi0 @ t00 @ Phi0 @ t0t @ Phit @ ttt @ Phit @ tt0)
        VS_t0 = np.trace(Phit @ ttt @ Phit @ tt0 @ Phi0 @ t0t)
        diags[0, dt] = -SS_0t
        diags[1, dt] = -TS
        diags[2, dt] = DSt * VS_0t
        diags[3, dt] = DS0 * VS_t0
        diags[4, dt] = CS ** 2
        diags[5, dt] = DpS0 * DpSt
        diags[6, dt] = -DS0 * DSt * CS
        diags[7, dt] = -DS0 ** 2 * DpSt
        diags[8, dt] = -DSt ** 2 * DpS0
        diags[9, dt] = DS0 ** 2 * DSt ** 2
    return diags


def build_C(Phi, tau, taugw, twin):
    # C[a,b,dt] for a,b in {PS,FS} = S_4 (all tau) + S~_4 (sink-a leg at t, sink-b leg at 0)
    legs = [tau, -taugw]
    S4_plain = np.tensordot(dc.W10, compute_diags2(Phi, tau, tau, twin), axes=(0, 0))
    C = np.zeros((2, 2, twin), complex)
    for a in range(2):
        for b in range(2):
            St = np.tensordot(dc.W10, compute_diags2(Phi, legs[a], legs[b], twin), axes=(0, 0))
            C[a, b] = S4_plain + St
    return C


def gevp_energies(Cmat, t0):
    # asymmetric GEVP C(t) v = lambda C(t0) v -> E_n = -log|lambda_n| / (t-t0), ground = largest |lambda|
    twin = Cmat.shape[2]
    C0 = Cmat[:, :, t0]
    n = C0.shape[0]
    E = np.full((twin, n), np.nan)
    for t in range(t0 + 1, twin):
        try:
            Minv = np.linalg.solve(C0, Cmat[:, :, t])
            lam = np.linalg.eigvals(Minv)
            lam = lam[np.argsort(-np.abs(lam))]         # largest |lambda| first = lowest E
            E[t] = -np.log(np.abs(lam)) / (t - t0)
        except Exception:
            pass
    return E


def main():
    print("# ENS=%s  ncfg=%d" % (dc.ENS.split("nu0")[0], len(dc.KS)))
    if len(dc.KS) < 10:
        print("# too few configs")
        return
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    Cs = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        Cs.append(build_C(Phi, tau, taugw, twin))
    Cs = np.array(Cs)                      # (ncfg,2,2,twin)
    ncfg = Cs.shape[0]
    twin = Cs.shape[3]
    Cavg = Cs.mean(0).real                 # (2,2,twin)

    # asymmetry check
    asy = np.abs(Cavg[0, 1] - Cavg[1, 0]) / (np.abs(Cavg[0, 1]) + np.abs(Cavg[1, 0]) + 1e-300)
    print("# PS-FS asymmetry |C01-C10|/|C01+C10| at dt=1,4,8 = %.3f %.3f %.3f" % (asy[1], asy[4], asy[8]))

    # vacuum subtraction: subtract the large-dt plateau (= <O_a><O_b>)
    plat = Cavg[:, :, twin - 6:].mean(2)
    Cc = Cavg - plat[:, :, None]

    # asymmetric GEVP at a few t0
    for t0 in (2, 3, 4):
        E = gevp_energies(Cc, t0)
        print("\n=== asymm GEVP, t0=%d  (E0=ground, E1=excited) ===" % t0)
        print("  t     E0        E1")
        for t in range(t0 + 1, min(20, twin)):
            print("  %2d   %7.4f   %7.4f" % (t, E[t, 0], E[t, 1]))
    # ground-state eigenvector content at a representative (t,t0)
    t0, t = 3, 10
    C0 = Cc[:, :, t0]
    vR = np.linalg.eig(np.linalg.solve(C0, Cc[:, :, t]))[1]
    v0 = vR[:, np.argmax(np.abs(np.linalg.eigvals(np.linalg.solve(C0, Cc[:, :, t]))))]
    v0 = v0 / np.abs(v0).max()
    print("\n# ground eigenvector (t=%d,t0=%d) content: PS=%.3f  FS=%.3f (|v| normalized)"
          % (t, t0, abs(v0[0]), abs(v0[1])))


def load_glue_OF(k, twin):
    # l=0 F^2 per-timeslice operator O_F(s) (glue op 0), restricted to the window [0,twin)
    import h5py
    gl = "data_" + dc.ENS + "/glue_f2_v2_shapes.%d.h5" % k
    with h5py.File(gl, "r") as h:
        O = np.array(h["O"])[0]                 # (Nt,)
    return O[:twin].astype(float)


def sigma_op_values(Phi, tau, taugw, twin):
    # per-config gauge-functional operator values O_a(s) = D_S,a(s)^2 + D'_S,a(s) for a in {PS,FS}.
    # sigma_a = eta^dag S xi + xi^dag S~ eta -> (S+S~) leg = tau + leg_S~ ; leg_S~: PS=tau, FS=-taugw.
    OPS = np.zeros(twin)
    OFS = np.zeros(twin)
    for s in range(twin):
        Ps = Phi[s]
        lsum_ps = 2.0 * tau[s, s]                       # (S+S~) for PS = 2 tau
        lsum_fs = tau[s, s] - taugw[s, s]               # (S+S~) for FS = tau - tau'
        DS_ps = np.trace(Ps @ lsum_ps).real
        DS_fs = np.trace(Ps @ lsum_fs).real
        Dp_ps = np.trace(Ps @ lsum_ps @ Ps @ lsum_ps).real
        Dp_fs = np.trace(Ps @ lsum_fs @ Ps @ lsum_fs).real
        OPS[s] = DS_ps ** 2 + Dp_ps
        OFS[s] = DS_fs ** 2 + Dp_fs
    return OPS, OFS


def full_gevp():
    # 3x3 mixing GEVP {F^2, sigma_PS^2, sigma_FS^2} from per-config operator values (gauge-functional sigma).
    print("\n########## CHUNK 2: F^2-sigma\\sigma mixing GEVP (3x3, gauge-functional sigma) ##########")
    print("# ENS=%s  ncfg=%d" % (dc.ENS.split("nu0")[0], len(dc.KS)))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    labels = ["F2", "PS2", "FS2"]
    Oall = []
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Phi = []
        for a in range(twin):
            Vt = V[tsrc0 + a].T
            Phi.append(Vt.conj().T @ (w00[:, None] * Vt))
        OF = load_glue_OF(k, twin)
        OPS, OFS = sigma_op_values(Phi, tau, taugw, twin)
        Oall.append(np.stack([OF, OPS, OFS]))          # (3, twin)
    Oall = np.array(Oall)                              # (ncfg, 3, twin)
    ncfg, nop, twin = Oall.shape

    def corr_matrix(O):
        # C[a,b,dt] = <O_a(s+dt) O_b(s)> - <O_a><O_b>, TRANSLATION-AVERAGED over window sources s.
        nc, no, tw = O.shape
        Obar = O.mean(axis=(0, 2))
        C = np.zeros((no, no, tw))
        for dt in range(tw):
            ns = tw - dt
            acc = np.zeros((nc, no, no))
            for s in range(ns):
                acc += O[:, :, s + dt][:, :, None] * O[:, :, s][:, None, :]
            acc /= ns
            C[:, :, dt] = acc.mean(0) - np.outer(Obar, Obar)
        return C

    def gevp_full(C, t0):
        twin = C.shape[2]
        C0 = C[:, :, t0]
        E = np.full((twin, C0.shape[0]), np.nan)
        vecs = {}
        for t in range(t0 + 1, twin):
            try:
                Minv = np.linalg.solve(C0, C[:, :, t])
                lam, vr = np.linalg.eig(Minv)
                order = np.argsort(-np.abs(lam))
                lam = lam[order]
                vr = vr[:, order]
                E[t] = -np.log(np.abs(lam)) / (t - t0)
                vecs[t] = vr
            except Exception:
                pass
        return E, vecs

    C = corr_matrix(Oall)
    # F^2 and sigma effmass sanity (diagonal)
    for a in range(nop):
        with np.errstate(all="ignore"):
            m = np.log(C[a, a, :-1] / C[a, a, 1:])
        print("  [%s diag] effmass dt=2,4,6,8 = %6.3f %6.3f %6.3f %6.3f" % (labels[a], m[2], m[4], m[6], m[8]))
    t0 = 3
    E, vecs = gevp_full(C, t0)
    print("\n  === 3x3 GEVP energies (t0=%d) ===" % t0)
    print("  t     E0       E1       E2")
    for t in range(t0 + 1, min(16, twin)):
        print("  %2d  %7.4f  %7.4f  %7.4f" % (t, E[t, 0], E[t, 1], E[t, 2]))
    # light-branch eigenvector content (jackknife-free central value)
    t = 8
    if t in vecs:
        v0 = np.abs(vecs[t][:, 0])
        v0 = v0 / v0.max()
        print("\n  # light-branch (E0) eigenvector content @t=%d: F2=%.3f PS2=%.3f FS2=%.3f" % (t, v0[0], v0[1], v0[2]))
        print("  # -> is the light state a MIX of F^2 and sigma\\sigma?  (Delta_- interpretation)")
    # jackknife on E0 plateau
    E0s = []
    for i in range(ncfg):
        Oj = np.delete(Oall, i, axis=0)
        Cj = corr_matrix(Oj)
        Ej, _ = gevp_full(Cj, t0)
        E0s.append(Ej[:, 0])
    E0s = np.array(E0s)
    E0m = E0s.mean(0)
    E0e = np.sqrt((ncfg - 1) * np.mean((E0s - E0m) ** 2, 0))
    print("\n  === E0 (ground) with jackknife errors, t0=%d ===" % t0)
    for t in range(t0 + 1, min(14, twin)):
        print("  t=%2d  E0 = %.4f(%.4f)" % (t, E0m[t], E0e[t]))


if __name__ == "__main__":
    main()
    full_gevp()
