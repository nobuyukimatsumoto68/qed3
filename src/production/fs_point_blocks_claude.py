#!/usr/bin/env python3
# fs_point_blocks_claude.py  [CHUNK 1: position-space point-block machinery + Y00 validation]
# Run:  ENS=... NVDIR=distill_Nv24 python3 fs_point_blocks_claude.py
#
# Position-space propagator  A(t,s) = U_t tau(t,s) U_s^dag   (2Ns x 2Ns ; EXACT D^-1 at L1 where U U^dag=1),
#   U_t = V[t].T  (2Ns x Nv, columns = distillation vectors).  Site block A(t,s)_{xy} = A[2x:2x+2, 2y:2y+2].
# Point building blocks (spin-traced 2x2 blocks):
#   point meson   G_xy(t,s) = -Tr_spin[ A(t,s)_{xy} A(s,t)_{yx} ]           = -Tr_Nv[Phi_x tau Phi_y tau]
#   point tadpole DS(x,t)   =  Tr_spin[ (A(t,t)-contact)_{xx} ]             =  Tr_Nv[Phi_x tt]
# CONTACT (equal time): S-part A_S(t,t) - 1/2 I ; Stilde  -1/2 (A'(t,t)+A(t,t)) ; A'(t,s)=-U_t tau'(t,s) U_s^dag.
# VALIDATION: the Y00 area-sum of point blocks must reproduce the Nv-basis Phi_00 quantities:
#   C_S^00(t,s) = sum_{x,y} A_x A_y Y00^2 G_xy(t,s)   ==  -Tr[Phi_00 tau Phi_00 tau]
#   DS^00       = sum_x A_x Y00 DS(x,t)               ==   Tr[Phi_00 tt]
# (leg tau = S-part; -tau' = Stilde-part.)  See fs_gevp_connected_impl_plan_claude.md.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

NS = dc.NS


def posprop(U, tau, t, s):
    # A(t,s) = U_t tau(t,s) U_s^dag   (2Ns x 2Ns)
    return U[t] @ tau[t, s] @ U[s].conj().T


def site_blocks(A):
    # reshape (2Ns,2Ns) -> (Nsite, 2, Nsite, 2)
    n2 = A.shape[0]
    ns = n2 // NS
    return A.reshape(ns, NS, ns, NS)


def main():
    tag = dc.ENS.split("nu0")[0]
    dual = dc.dual_areas_from_mesh()
    nsite = dual.shape[0]
    w00 = np.repeat(dual, NS) * dc.Y00
    k = dc.KS[0]
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Nv = V.shape[1]
    Iv = np.eye(Nv)
    # U_t = V[t].T  (2Ns x Nv), aligned to the window index a -> timeslice tsrc0+a
    U = [V[tsrc0 + a].T for a in range(twin)]
    # Nv-basis reference vertices
    Phi = [(U[a].conj().T) @ (w00[:, None] * U[a]) for a in range(twin)]

    # completeness check (L1: U U^dag = I_{2Ns})
    UU = U[0] @ U[0].conj().T
    print("# ENS=%s  Nv=%d 2Ns=%d nsite=%d  ||U U^dag - I|| = %.2e (0 => complete/exact posprop)"
          % (tag, Nv, U[0].shape[0], nsite, np.linalg.norm(UU - np.eye(UU.shape[0]))))

    t, s = 6, 0
    # ---- S-part (leg tau) ----
    A_ts = posprop(U, tau, t, s)
    A_st = posprop(U, tau, s, t)
    Bts = site_blocks(A_ts)
    Bst = site_blocks(A_st)
    # point meson G_xy(t,s) = -Tr_spin[ A(t,s)_xy A(s,t)_yx ]
    G = -np.einsum("xayb,ybxa->xy", Bts, Bst).real          # (nsite,nsite)
    CS00_point = np.einsum("x,y,xy->", dual * dc.Y00, dual * dc.Y00, G)
    CS00_nv = (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
    print("# [S-part] Y00 meson C_S(t=%d,s=%d): point-sum=% .8e  Nv=% .8e  ratio=%.8f"
          % (t, s, CS00_point, CS00_nv, CS00_point / CS00_nv))

    # point tadpole DS(x,t) = Tr_spin[ (A(t,t)-1/2 I)_xx ]
    A_tt = posprop(U, tau, t, t) - 0.5 * np.eye(A_ts.shape[0])
    Btt = site_blocks(A_tt)
    DS = np.einsum("xaxa->x", Btt).real
    DS00_point = np.einsum("x,x->", dual * dc.Y00, DS)
    tt_nv = tau[t, t] - 0.5 * Iv
    DS00_nv = np.trace(Phi[t] @ tt_nv).real
    print("# [S-part] Y00 tadpole D_S(t=%d): point-sum=% .8e  Nv=% .8e  (both ~0 = contact)"
          % (t, DS00_point, DS00_nv))

    # point self-loop D'_S = sum_{x,x'} (A_x Y00)(A_x' Y00) Tr_spin[ tt_{x x'} tt_{x' x} ]  (double site sum)
    w = dual * dc.Y00
    DpS00_point = np.einsum("x,y,xayb,ybxa->", w, w, Btt, Btt).real
    DpS00_nv = np.trace(Phi[t] @ tt_nv @ Phi[t] @ tt_nv).real
    print("# [S-part] Y00 self-loop D'_S(t=%d): point-sum=% .8e  Nv=% .8e  ratio=%.8f"
          % (t, DpS00_point, DpS00_nv, DpS00_point / DpS00_nv))

    # ---- Stilde-part (leg -tau') : FS contact -1/2 (A'(t,t)+A(t,t)) ----
    Ap_ts = -posprop(U, taugw, t, s)
    Ap_st = -posprop(U, taugw, s, t)
    Gp = -np.einsum("xayb,ybxa->xy", site_blocks(Ap_ts), site_blocks(Ap_st)).real
    CSp_point = np.einsum("x,y,xy->", dual * dc.Y00, dual * dc.Y00, Gp)
    # FS leg collapses to tau by GW (S~ D_ov^{-dag} = tau); -taugw was the tau_gw artifact -> FS == PS.
    # See fs_furnishing_derivation_claude.md.  Original (A/B):
    # CSp_nv = (-np.trace(Phi[t] @ (-taugw[t, s]) @ Phi[s] @ (-taugw[s, t]))).real
    CSp_nv = (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
    print("# [Stilde] Y00 meson C_S[-tau'](t=%d,s=%d): point-sum=% .8e  Nv=% .8e  ratio=%.8f"
          % (t, s, CSp_point, CSp_nv, CSp_point / CSp_nv))


if __name__ == "__main__":
    main()
