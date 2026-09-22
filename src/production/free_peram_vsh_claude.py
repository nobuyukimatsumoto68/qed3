# Validate the peram->VSH pipeline: reconstruct the (free) all-to-all overlap propagator from the
# distillation perambulator (complete basis), build the tangential current tensor f^{ab}(n1,n2), project
# onto VSH (electric Phi / magnetic Psi, ell=1), and compare to the direct prop_deter result
# (free_L1_vsh_claude.py: L1 lattice E~0.38, M~0.42; continuum 0.40/0.60).
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ["ENS"] = "free"
os.environ["NVDIR"] = "distill_Nv24"
import sys
import numpy as np
import fs_gevp_point_claude as fg
import distill_contract_claude as dc
sys.path.insert(0, "final/analysis_axial")
import effmass_axial_tp_l3_perm_hankel_claude as H

GEO = "/mnt/barracuda22/qed3/qed3/geometry/data"
AT = 0.2
L = 1
s1 = np.array([[0, 1], [1, 0]], complex)
s2 = np.array([[0, -1j], [1j, 0]], complex)
SIG = {1: s1, 2: s2}


def ylm1_weights(theta, phi):
    c = np.sqrt(3.0 / (4.0 * np.pi))
    st = np.sin(theta)
    ct = np.cos(theta)
    return {0: (-c * st, np.zeros_like(st)),
            1: (c * ct * np.cos(phi), -c * np.sin(phi)),
            -1: (c * ct * np.sin(phi), c * np.cos(phi))}


def main():
    AblkS, AblkSt, twin, nsite, U, tau = fg.make_config(0)
    pts = np.loadtxt("%s/pts_n%d.dat" % (GEO, L))[:nsite]
    n = pts / np.linalg.norm(pts, axis=1, keepdims=True)
    theta = np.arccos(np.clip(n[:, 2], -1, 1))
    phi = np.arctan2(n[:, 1], n[:, 0])
    w = np.ones(nsite)
    Wt = ylm1_weights(theta, phi)
    print("# free peram: twin=%d nsite=%d  (VSH ell=1; targets lattice E~0.38 M~0.42, continuum 0.40/0.60)" % (twin, nsite))

    dmax = twin - 1
    gE = np.zeros(dmax)
    gM = np.zeros(dmax)
    for dt in range(1, dmax + 1):
        eE = eM = 0.0
        ncnt = 0
        for t0 in range(0, twin - dt):
            G12 = AblkS(t0, t0 + dt).transpose(0, 2, 1, 3)          # [n1,n2,a,b] = S(n1@t0, n2@t0+dt)
            G21 = AblkS(t0 + dt, t0).transpose(2, 0, 1, 3)          # [n1,n2,a,b] = S(n2@t0+dt, n1@t0)
            f = {}
            for a in (1, 2):
                for b in (1, 2):
                    f[(a, b)] = -np.einsum('pq,ijqr,rs,ijsp->ij', SIG[a], G12, SIG[b], G21, optimize=True)
            for m in (-1, 0, 1):
                wt, wp = Wt[m]
                WE = {1: wt, 2: wp}
                WM = {1: -wp, 2: wt}
                for a in (1, 2):
                    for b in (1, 2):
                        eE += np.sum(np.outer(w * WE[a], w * WE[b]) * f[(a, b)])
                        eM += np.sum(np.outer(w * WM[a], w * WM[b]) * f[(a, b)])
            ncnt += 1
        gE[dt - 1] = (eE / ncnt).real
        gM[dt - 1] = (eM / ncnt).real

    off = [0, 3, 6]
    for lab, C in [("E (electric)", gE), ("M (magnetic)", gM)]:
        sgn = np.sign(C[3])
        emH, _ = H.hankel_effmass_scalar((sgn * C).astype(float), off, 4, 1, 1, AT)
        naive = np.log(np.abs(C[:-1] / C[1:]))
        print(" %-13s GPOF =%s" % (lab, np.array2string(emH[:14, 0], precision=3, floatmode="maxprec")))
        print("               naive=%s" % np.array2string(naive[:14], precision=3, floatmode="maxprec"))


if __name__ == "__main__":
    main()
