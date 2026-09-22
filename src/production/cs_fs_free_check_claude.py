#!/usr/bin/env python3
# cs_fs_free_check_claude.py
#   Free-limit sign/orientation check of taugw (Fin's validation (2), free version): the FS single-meson
#   two-point built from the furnished leg -taugw must plateau at the single-meson mass 2E0 (=0.378 L1).
#   C_S^FS(dt) = < -Tr[ Phi_t (-taugw[t,s]) Phi_s (-taugw[s,t]) ] >_s .  If effmass -> 0.378, taugw is right.
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 cs_fs_free_check_claude.py

import os
os.environ.setdefault("ENS", "free")
os.environ.setdefault("LREF", "1")
os.environ.setdefault("NVDIR", "distill_Nv24")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc


def main():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    C = np.full(twin, np.nan)
    for dt in range(1, twin):
        ns = twin - dt
        if ns <= 0:
            continue
        acc = 0.0
        for s in range(ns):
            t = s + dt
            acc += (-np.trace(Phi[t] @ (-taugw[t, s]) @ Phi[s] @ (-taugw[s, t]))).real
        C[dt] = acc / ns

    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    with np.errstate(all="ignore"):
        em = np.log(np.abs(C[:-1] / C[1:]))
    print("# FREE L=%d  C_S^FS (from -taugw) effmass -- should plateau at 2E0 = %.3f (taugw sign check)" % (dc.L, msig))
    print("#  dt |   C_S^FS         effmass")
    for dt in range(1, twin - 1):
        if np.isfinite(C[dt]):
            e = em[dt] if np.isfinite(em[dt]) else np.nan
            print("#  %2d | %+.6e   %7.4f" % (dt, C[dt], e))


if __name__ == "__main__":
    main()
