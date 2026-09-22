#!/usr/bin/env python3
# triangle_split_claude.py -- WHY does the cross C_12 = <sigma_00(t) sigma^2_00(0)> vanish?  Per-loop vanishing,
#   or cancellation between the two loop orientations?  The triangle has two orientations (the two source vertices
#   k,l swapped): cyc1 = [0,1,2] (i->k->l), cyc2 = [0,2,1] (i->l->k).  This splits C_12 into T1 (cyc1 only),
#   T2 (cyc2 only), and reports each ABSOLUTE (not normalized) plus the sum, at several dt, for whatever ENS/NVDIR
#   is set.  Also the per-vertex-site DIAGONAL vs full (is the spinor trace zero site-by-site, or only after the
#   area sum?).  Reuses sigma2_mPS_gevp.perm_contrib_n / G.make_config_win / G.op_vspec.
#   Run: ENS=.. LREF=.. NVDIR=.. [NVKEEP=..] [MODE_CONTACT=..] python3 triangle_split_claude.py

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
DTS = [int(x) for x in os.environ.get("DTS", "2,4,6,8").split(",")]


def main():
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = G.antipodal_map()
    AblkS, AblkSt, twin, nsite, U, tau, tsrc0 = G.make_config_win(KCFG, 0)
    vsp_12 = [('i', wY, False)] + G.op_vspec(0, ('k', 'l'), dualf, wY)   # sink sigma_00 (i) + source sigma^2 (k,l)
    vsp_ss = [(('i', wY, False)), (('k', wY, False))]

    print("# ENS=%s L=%d NVDIR=%s NVKEEP=%s MODE_CONTACT=%s k=%d"
          % (dc.ENS.split("nu0")[0], dc.L, os.environ["NVDIR"], os.environ.get("NVKEEP", "0"),
             os.environ.get("MODE_CONTACT", "0"), KCFG))
    print("# triangle C_12 = T1 + T2 :  T1 = cyc[0,1,2] only,  T2 = cyc[0,2,1] only  (each = single orientation)")
    print("# dt |        T1              T2           T1+T2 (=C12)      C11        |T1+T2|/|C11|  |T1-T2|/(|T1|+|T2|)")
    for dt in DTS:
        s0s = np.array([s for s in range(twin) if s + dt < twin])
        offs = {(0, 0), (dt, dt), (0, dt), (dt, 0)}
        bAS = {o: np.array([AblkS(s + o[0], s + o[1]) for s in s0s]) for o in offs}
        T1 = MG.perm_contrib_n([[0, 1, 2]], [dt, 0, 0], bAS, vsp_12, Pmap).real / len(s0s)
        T2 = MG.perm_contrib_n([[0, 2, 1]], [dt, 0, 0], bAS, vsp_12, Pmap).real / len(s0s)
        C11 = MG.perm_contrib_n([[0, 1]], [dt, 0], bAS, vsp_ss, Pmap).real / len(s0s)
        s = T1 + T2
        rC = abs(s) / abs(C11) if C11 != 0 else np.nan
        rD = abs(T1 - T2) / (abs(T1) + abs(T2)) if (abs(T1) + abs(T2)) > 0 else np.nan
        print("# %2d | %14.6e  %14.6e  %14.6e  %12.4e  %11.3e   %11.3e"
              % (dt, T1, T2, s, C11, rC, rD))

    # site-by-site: is the SPINOR trace zero per (x,y,z), or only after the area sum?  Sample the raw loop tensor
    # for cyc[0,1,2] BEFORE the wY vertex weights, at dt = DTS[-1], and report ||diag(x=y=z)|| vs full mean.
    dt = DTS[-1]
    s0s = np.array([s for s in range(twin) if s + dt < twin])
    A_t0 = np.mean([AblkS(s + dt, s + 0) for s in s0s], axis=0)     # (nsite,NS,nsite,NS)
    A_00 = np.mean([AblkS(s + 0, s + 0) for s in s0s], axis=0)
    A_0t = np.mean([AblkS(s + 0, s + dt) for s in s0s], axis=0)
    # loop tensor L[x,y,z] = Tr_spin[ A_t0[x,:,y,:] A_00[y,:,z,:] A_0t[z,:,x,:] ]  (NO wY weights)
    L = np.einsum('xayb,ybzc,zcxa->xyz', A_t0, A_00, A_0t, optimize='optimal')
    diag = np.array([L[x, x, x] for x in range(nsite)])
    print("# [site-by-site, dt=%d, cyc[0,1,2], NO area weights] loop tensor L[x,y,z]=Tr_spin[...]:" % dt)
    print("#   max|L| (all x,y,z) = %.4e   mean|L| = %.4e   max|L[x,x,x]| (site-diag) = %.4e"
          % (np.max(np.abs(L)), np.mean(np.abs(L)), np.max(np.abs(diag))))
    print("#   -> if max|L| >> |C12| but the wY-weighted SUM ~0, the zero is a SUM cancellation, not site-by-site.")


if __name__ == "__main__":
    main()
