#!/usr/bin/env python3
# sigma_ps_fs_2pt_claude.py
#   Single-meson two-points: <sigma_PS sigma_PS>, <sigma_FS sigma_FS>, and the CROSS <sigma_PS sigma_FS>.
#   NM: is the PS<->FS cross zero?  Furnishing: leg arriving at an FS vertex -> -taugw; at PS -> plain tau.
#     C_PP = -Tr[Phi_t tau[t,s] Phi_s tau[s,t]]
#     C_FF = -Tr[Phi_t (-taugw[t,s]) Phi_s (-taugw[s,t])]
#     C_PF = -Tr[Phi_t tau[t,s] Phi_s (-taugw[s,t])]   (sink PS plain, source FS furnished)
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 sigma_ps_fs_2pt_claude.py

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

LTAG = os.environ.get("LTAG", "L%d" % dc.L)


def main():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    CPP = np.full(twin, np.nan)
    CFF = np.full(twin, np.nan)
    CPF = np.full(twin, np.nan)
    for dt in range(twin):
        ns = twin - dt
        if ns <= 0:
            continue
        pp = ff = pf = 0.0
        for s in range(ns):
            t = s + dt
            pp += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            ff += (-np.trace(Phi[t] @ (-taugw[t, s]) @ Phi[s] @ (-taugw[s, t]))).real
            pf += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ (-taugw[s, t]))).real
        CPP[dt], CFF[dt], CPF[dt] = pp / ns, ff / ns, pf / ns

    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# %s single-meson two-points; m_sig=%.3f.  R_PF = <PS FS>/<PS PS> (0 => cross vanishes)" % (LTAG, msig))
    print("#  dt |   C_PP          C_FF          C_PF(cross)    R_PF=C_PF/C_PP")
    for dt in range(1, twin - 1):
        if np.isfinite(CPP[dt]) and abs(CPP[dt]) > 0:
            print("#  %2d | %+.5e   %+.5e   %+.5e   %+.6f" % (dt, CPP[dt], CFF[dt], CPF[dt], CPF[dt] / CPP[dt]))


if __name__ == "__main__":
    main()
