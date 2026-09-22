#!/usr/bin/env python3
# tadpole_trunc_check_claude.py -- does BASIS TRUNCATION revive the contact-subtracted single-sigma tadpole
#   D_S(t) = Tr[Phi (tau(t,t) - 1/2 I)] ?  At the COMPLETE basis V V^dag = I so the GW contact 1/2 is exact and
#   D_S = 0 (config-indep).  Under truncation A(t,t) = U tau U^dag - 1/2 I with U U^dag = P != I -- the propagator
#   part lives in the truncated subspace but the 1/2 I contact is on the FULL 2Ns space, so the subtraction may be
#   inconsistent and D_S may go nonzero, feeding a single-meson piece into C_22 = <sigma^2 sigma^2> via D_S(t)D_S(0).
#   (If D_S stays ~0 under truncation, the C_22 collapse is NOT tadpole revival -> a connected/degenerate-state effect.)
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 NVLIST=0,18,12,6 python3 tadpole_trunc_check_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "free")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc

NVLIST = [int(x) for x in os.environ.get("NVLIST", "0,18,12,6").split(",")]
KCFG = int(os.environ.get("KCFG", str(dc.KS[0])))
NS = dc.NS


def main():
    dualf = dc.dual_areas_from_mesh().astype(float)
    w00 = np.repeat(dualf, NS) * dc.Y00        # per-spinor vertex weight (j = NS*x + s)
    print("# ENS=%s L=%d NVDIR=%s k=%d   contact tadpole D_S(t)=Tr[Phi(tau_tt - 1/2 I)] vs NVKEEP (free)"
          % (dc.ENS.split("nu0")[0], dc.L, os.environ["NVDIR"], KCFG))
    print("# Nv   |   mean_t D_S      max_t|D_S|     (0 = contact exact / no tadpole revival)")
    for nv in NVLIST:
        os.environ["NVKEEP"] = str(nv)
        V, windows = dc.load_peram_windows(KCFG)
        tsrc0, tau, taugw = windows[0]
        twin = tau.shape[0]
        twoNs = V.shape[2]
        Iv = np.eye(twoNs)
        nvk = tau.shape[-1]
        Inv = np.eye(nvk)
        DS = np.zeros(twin)                      # current code: contact -1/2 I on FULL 2Ns space
        DSm = np.zeros(twin)                     # FIX: contact -1/2 in MODE space (= -1/2 P), tau - 1/2 I_Nv
        for t in range(twin):
            Ut = V[tsrc0 + t].T                 # (2Ns, Nv)
            A = Ut @ tau[t, t] @ Ut.conj().T - 0.5 * Iv
            Am = Ut @ (tau[t, t] - 0.5 * Inv) @ Ut.conj().T
            DS[t] = np.real(np.sum(w00 * np.diag(A)))
            DSm[t] = np.real(np.sum(w00 * np.diag(Am)))
        lab = "all" if nv == 0 else "%d" % nv
        print("# %4s | %14.6e  %14.6e | %14.6e  %14.6e"
              % (lab, DS.mean(), np.max(np.abs(DS)), DSm.mean(), np.max(np.abs(DSm))))


if __name__ == "__main__":
    main()
