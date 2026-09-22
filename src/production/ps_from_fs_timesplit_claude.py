#!/usr/bin/env python3
# ps_from_fs_timesplit_claude.py
#   Fin's DECISIVE test: is the O(1) sigma_FS^2 -> single-meson coupling the FURNISHING (genuine) or the
#   unsolved EQUAL-TIME FS source contact (artifact)?  Time-split the FS source: src1 at t=0, src2 at t=Delta.
#   The src-src leg becomes a time-separated furnished propagator -taugw(0,Delta) with NO equal-time contact.
#   - persists smoothly for Delta=1,2 (extrapolating to Delta=0) -> FURNISHING structure, genuine.
#   - spikes at Delta=0, dies for Delta>0 -> equal-time contact was masquerading, artifact.
#   Channel: PS<-FS (sink sigma_PS, source sigma_FS^2).  Legs: arrival-at-FS -> -taugw; arrival-at-PS-sink -> tau.
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 ps_from_fs_timesplit_claude.py

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
DELTAS = [0, 1, 2]


def build():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    # sigma_PS sigma_PS reference two-point (for R normalization)
    Css = np.full(twin, np.nan)
    for dt in range(twin):
        ns = twin - dt
        if ns <= 0:
            continue
        acc = 0.0
        for s in range(ns):
            acc += (-np.trace(Phi[s + dt] @ tau[s + dt, s] @ Phi[s] @ tau[s, s + dt])).real
        Css[dt] = acc / ns

    # time-split PS<-FS triangle: src1 @ s (FS), src2 @ s+D (FS), sink @ s+dt (PS)
    #   loop  s -(L1)-> s+D -(L2)-> s+dt -(L3)-> s ;  arrival-FS -> -taugw, arrival-PS -> tau
    C = {D: np.full(twin, np.nan) for D in DELTAS}
    for D in DELTAS:
        for dt in range(D + 1, twin):
            ns = twin - dt
            if ns <= 0:
                continue
            acc = 0.0
            for s in range(ns):
                t = s + dt
                sd = s + D
                L1 = -taugw[sd, s]        # arrives src2 (FS)
                L2 = tau[t, sd]           # arrives sink (PS)
                L3 = -taugw[s, t]         # arrives src1 (FS)
                acc += (-np.trace(Phi[s] @ L1 @ Phi[sd] @ L2 @ Phi[t] @ L3)).real
            C[D][dt] = acc / ns
    return Css, C, twin


def main():
    Css, C, twin = build()
    msig = {1: 0.378, 2: 0.393}.get(dc.L, 0.0)
    print("# %s PS<-FS time-split: src1@0, src2@Delta.  R=C/Css (relative overlap), effmass; m_sig=%.3f" % (LTAG, msig))
    print("# genuine FURNISHING -> R persists for Delta>0; equal-time CONTACT artifact -> R dies for Delta>0.")
    for D in DELTAS:
        with np.errstate(all="ignore"):
            em = np.log(np.abs(C[D][:-1] / C[D][1:]))
        print("# --- Delta=%d ---" % D)
        for dt in range(max(D + 1, 6), min(twin - 2, 16)):
            if np.isfinite(Css[dt]) and abs(Css[dt]) > 0 and np.isfinite(C[D][dt]):
                e = em[dt] if np.isfinite(em[dt]) else np.nan
                print("#   dt=%2d  C=%+.4e  R=%+.5f  effmass=%7.4f" % (dt, C[D][dt], C[D][dt] / Css[dt], e))


if __name__ == "__main__":
    main()
