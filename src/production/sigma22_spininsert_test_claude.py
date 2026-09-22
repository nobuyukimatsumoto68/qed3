#!/usr/bin/env python3
# sigma22_spininsert_test_claude.py
#   Does inserting a spin matrix into the shell-projected (2,2) operator make it couple to sigma^2 (PS^2)?
#   O_22 = psibar Q2 psi is sigma3-EVEN -> <O_22 sigma^2>=0.  Test O_sp = psibar (sp) Q2 psi for
#   sp in {sigma1,sigma2,sigma3}: measure the connected mixing C2Q = <O_sp sigma^2> vs C2A = <O_A sigma^2>.
#   Coupling turns on only if the vertex is sigma3-ODD (the operative symmetry is sigma3 x coord-reflection,
#   so local spin insertions may stay even -- that is the point of the test).
#   Run:  ENS=free LREF=1 NVDIR=distill_Nv24 python3 sigma22_spininsert_test_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import diag_effmass_claude as de

DT_DECOMP = int(os.environ.get("DT_DECOMP", "3"))
T0 = int(os.environ.get("T0", "3"))


def shell_projector(tau, twin, Nv, dt_dec):
    ns = twin - dt_dec
    K = np.zeros((Nv, Nv), complex)
    for s in range(ns):
        K += tau[s + dt_dec, s]
    K /= ns
    mu, R = np.linalg.eig(K)
    Rinv = np.linalg.inv(R)
    E = -np.log(np.abs(mu)) / dt_dec
    order = np.argsort(E)
    R = R[:, order]
    Rinv = Rinv[order, :]
    E = E[order]
    clusters = []
    i = 0
    while i < Nv:
        j = i
        while j + 1 < Nv and abs(E[j + 1] - E[i]) < 0.02:
            j += 1
        clusters.append((i, j))
        i = j + 1
    i2, j2 = clusters[1]
    sel = list(range(i2, j2 + 1))
    return R[:, sel] @ Rinv[sel, :]


def main():
    de.CONTACT = 0.5
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s FREE L=%d  spin-insert coupling test for shell-projected (2,2)" % (tag, dc.L))
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Nv = V.shape[1]
    n2 = V.shape[2]
    nsite = n2 // dc.NS
    Iv = np.eye(Nv)
    tt = tau.copy()
    for a in range(twin):
        tt[a, a] = tau[a, a] - 0.5 * Iv
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]
    PA = [Phi[a] @ tt[a, a] for a in range(twin)]
    Q2 = shell_projector(tau, twin, Nv, DT_DECOMP)

    # spin matrices in the 2Ns space (index j = NS*x + s), block-diagonal over sites
    s0 = np.array([[1, 0], [0, 1]], complex)
    s1 = np.array([[0, 1], [1, 0]], complex)
    s2 = np.array([[0, -1j], [1j, 0]], complex)
    s3 = np.array([[1, 0], [0, -1]], complex)
    def spin_mode(sp, t):
        Sfull = np.kron(np.eye(nsite), sp)                 # (2Ns,2Ns)
        Vt = V[tsrc0 + t].T                                # (2Ns,Nv)
        return Vt.conj().T @ Sfull @ Vt                    # (Nv,Nv)

    # reference: C2A (O_A mixing with sigma^2) at a few t
    def C2_mix(Gam_list):
        out = np.zeros(twin)
        for dt in range(twin):
            ns = twin - dt
            acc = 0.0
            for s in range(ns):
                t = s + dt
                acc += (-2.0 * np.trace(Phi[s] @ tt[s, s] @ Phi[s] @ tau[s, t] @ Gam_list[t] @ tau[t, s])).real
            out[dt] = acc / ns
        return out
    def onepoint(Gam_list):
        return np.mean([(-np.trace(Gam_list[a] @ tt[a, a])).real for a in range(twin)])

    C2A = C2_mix(PA)
    oA = onepoint(PA)
    print("# reference O_A:  <O_A>=%.4e   |C2A(t0=%d)|=%.4e" % (oA, T0, abs(C2A[T0])))

    print("#\n#  insertion (order)     <O>          |C2Q(t0)|      |C2Q/C2A|")
    for name, sp in [("s3", s3), ("s1", s1), ("s2", s2)]:
        for order in ("sp@Q2", "Q2@sp"):
            Gam = []
            for t in range(twin):
                Sm = spin_mode(sp, t)
                Gam.append(Sm @ Q2 if order == "sp@Q2" else Q2 @ Sm)
            C2Q = C2_mix(Gam)
            oQ = onepoint(Gam)
            print("#  %s  %-8s      %+.3e   %.4e   %.3e"
                  % (name, order, oQ, abs(C2Q[T0]), abs(C2Q[T0]) / (abs(C2A[T0]) + 1e-300)))


if __name__ == "__main__":
    main()
