#!/usr/bin/env python3
# t00_cross_dump_claude.py
# Dump the free-L1 T_00 SINK vertex elementals Phi_W(a), Phi_WH(a) for the {1,1,1,1} thread's
# cross-correlator <T_00(t) sigma^2_P+(0)>.  Same peram window + V basis that distill_contract_claude.load_peram
# uses (ENS=free, NVDIR=distill_Nv24), so the sink shares the eigenvector basis and time index with the sigma^2 source.
#   O_H(t) = eta^H W(t) xi + xi^H W(t)^H eta  (r=0 naive e.sigma spatial hop; TWO terms)
#   Phi_W(a)  = V[tsrc0+a]^dag W  V[tsrc0+a]      (pairs with tau legs)
#   Phi_WH(a) = V[tsrc0+a]^dag W^dag V[tsrc0+a]   (h.c. term)
# Free: spatial link phase theta=0, so W = the ungauged geometry hop (identical to t00_stress_ham build).
# Run:  OMP_NUM_THREADS=4 python3 t00_cross_dump_claude.py

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ["ENS"] = "free"
os.environ["NVDIR"] = "distill_Nv24"
os.environ.setdefault("LREF", "1")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc
import geom_hopping_claude as gh
import t00_stress_ham_interacting_claude as ti

LREF = int(os.environ.get("LREF", "1"))
GEOM = os.environ.get("GEOM", "../../geometry/data/")
OUTDIR = os.environ.get("OUTDIR", "t00_cross_claude")


def main():
    om, alpha, nns, nsite = gh.build(GEOM, LREF)
    tab = ti.load_link_table("primal_links_n%d_claude.dat" % LREF)
    n_links = len(tab) // 2
    k = dc.KS[0]
    V, tau, taugw, tsrc0, twin = dc.load_peram(k)
    Nv = V.shape[1]
    Nt = V.shape[0]
    sp0 = np.zeros((Nt, n_links))                      # free: theta_sp = 0 -> phase = 1
    phiW = np.zeros((twin, Nv, Nv), complex)
    phiWH = np.zeros((twin, Nv, Nv), complex)
    for a in range(twin):
        t = tsrc0 + a
        Vt = V[t].T                                    # (2Ns, Nv)
        W = ti.build_W_gauge(om, alpha, nns, nsite, tab, sp0[t])   # r=0 by module default
        phiW[a] = Vt.conj().T @ W @ Vt
        phiWH[a] = Vt.conj().T @ W.conj().T @ Vt
    os.makedirs(OUTDIR, exist_ok=True)
    out = "%s/phiW_free_L1_claude.npz" % OUTDIR
    np.savez(out, phiW=phiW, phiWH=phiWH)
    print("# wrote %s" % out)
    print("# keys: phiW, phiWH   shape=(twin,Nv,Nv)=(%d,%d,%d)" % (twin, Nv, Nv))
    print("# ENS=free NVDIR=distill_Nv24 k=%d tsrc0=%d twin=%d Nv=%d nsite=%d n_links=%d r=%.1f"
          % (k, tsrc0, twin, Nv, nsite, n_links, ti.R))
    print("# index: a=0..twin-1 maps to timeslice tsrc0+a; Phi_W(a)=V[tsrc0+a]^dag W V[tsrc0+a]")
    print("# leg convention: <xi eta^dag>=tau, <eta xi^dag>=delta-tau; both sink terms tie in with tau")


if __name__ == "__main__":
    main()
