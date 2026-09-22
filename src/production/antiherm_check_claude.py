#!/usr/bin/env python3
# antiherm_check_claude.py -- confirm the mechanism NM proposed: the cross triangle vanishes config-by-config as
#   a loop identity Re Tr[(D_ov^{-1} - 1/2)^3] = 0, because the NORMAL-ORDERED propagator M = D_ov^{-1} - 1/2 delta
#   is ANTI-HERMITIAN (GW: D_ov^{-1} + D_ov^{-dag} = 1 => M^dag = D_ov^{-dag} - 1/2 = 1/2 - D_ov^{-1} = -M), so any
#   ODD closed loop of M (3 vertices) has T^* = -T => purely imaginary => Re = 0.  Checks, config-by-config:
#     (a) equal-time tt = tau(0,0) - 1/2 I_Nv is anti-hermitian:  ||tt + tt^dag|| / ||tt||  (should be ~0).
#     (b) off-diagonal-time GW:  tau(a,b)^dag = -tau(b,a)  (a!=b):  ||tau_ab^dag + tau_ba|| / ||tau_ab||.
#     (c) the triangle T1 is PURELY IMAGINARY: |Re T1| / |Im T1|  (should be ~0), while Im T1 is O(1)-ish.
#   Run: ENS=.. LREF=.. NVDIR=.. python3 antiherm_check_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import sigma2_mPS_gevp_claude as MG

G = MG.G
KCFG = int(os.environ.get("KCFG", str(dc.KS[0])))
DTS = [int(x) for x in os.environ.get("DTS", "2,4,6").split(",")]


def main():
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    V, tau, taugw, tsrc0, twin = dc.load_peram(KCFG)
    nvk = tau.shape[-1]
    Inv = np.eye(nvk)
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s L=%d NVDIR=%s k=%d  Nv=%d  -- anti-hermiticity of M = D_ov^{-1} - 1/2 (GW), config-by-config"
          % (tag, dc.L, os.environ["NVDIR"], KCFG, nvk))

    # (a) equal-time contact-subtracted block, mode space
    tt00 = tau[0, 0] - 0.5 * Inv
    ah = np.linalg.norm(tt00 + tt00.conj().T) / np.linalg.norm(tt00)
    print("# (a) equal-time tt(0,0)=tau-1/2 anti-herm:  ||tt+tt^dag||/||tt|| = %.3e   (0 => anti-hermitian)" % ah)

    # (b) off-diagonal-time GW:  tau(a,b)^dag = -tau(b,a)
    for (a, b) in [(2, 0), (4, 0), (6, 2)]:
        if a < twin and b < twin:
            r = np.linalg.norm(tau[a, b].conj().T + tau[b, a]) / np.linalg.norm(tau[a, b])
            print("# (b) off-diag tau(%d,%d)^dag + tau(%d,%d):  ||.||/||tau|| = %.3e   (0 => GW M^dag=-M)"
                  % (a, b, b, a, r))

    # (c) the triangle T1 (cyc[0,1,2]) real vs imag, per dt
    AblkS, AblkSt, twin2, nsite, U, tau2, tsrc0b = G.make_config_win(KCFG, 0)
    vsp_12 = [('i', wY, False)] + G.op_vspec(0, ('k', 'l'), dualf, wY)
    print("# (c) triangle T1 (cyc[0,1,2]) -- should be PURELY IMAGINARY (Re=0):")
    print("#  dt |     Re T1          Im T1        |Re|/|Im|")
    for dt in DTS:
        s0s = np.array([s for s in range(twin2) if s + dt < twin2])
        offs = {(0, 0), (dt, dt), (0, dt), (dt, 0)}
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        T1 = MG.perm_contrib_n([[0, 1, 2]], [dt, 0, 0], bAS, vsp_12, Pmap) / len(s0s)
        ratio = abs(T1.real) / abs(T1.imag) if T1.imag != 0 else np.nan
        print("#  %2d | %14.6e  %14.6e   %.3e" % (dt, T1.real, T1.imag, ratio))


if __name__ == "__main__":
    main()
