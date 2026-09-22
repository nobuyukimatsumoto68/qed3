#!/usr/bin/env python3
# flavor_triangle_duality_claude.py
#   All four single-meson triangles <sigma_X(t) sigma_Y^2(0)>, X,Y in {PS,FS}, single free config.
#   Furnishing rule (Fin, fs_ps2_triangle_recipe): each leg is furnished (-> -taugw) iff its ARRIVAL (row)
#   vertex is FS; the equal-time source self-leg uses the PS contact tau-1/2 (Y=PS) or the FS contact
#   -1/2(taugw+tau) (Y=FS).  Loop: src1 -(tt)-> src2 -(src_out, arrives src2)-> ... actually the three legs
#   arrive at src1(tt), sink(sink_in), src2(src_out); furnish each by its arrival flavor.
#   Expected duality: PS->PS = 0 (sigma3-herm) and FS->FS = 0 (FS-loop dual); PS->FS and FS->PS nonzero (O(a)).
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 flavor_triangle_duality_claude.py

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


def build():
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual.astype(float), dc.NS) * dc.Y00
    V, tau, taugw, tsrc0, twin = dc.load_peram(dc.KS[0])
    Iv = np.eye(tau.shape[-1])
    Phi = [(V[tsrc0 + a].T).conj().T @ (w00[:, None] * (V[tsrc0 + a].T)) for a in range(twin)]

    chans = [("PS", "PS"), ("PS", "FS"), ("FS", "PS"), ("FS", "FS")]   # (sink X, source Y)
    C = {c: np.full(twin, np.nan) for c in chans}
    Css = np.full(twin, np.nan)
    for dt in range(twin):
        ns = twin - dt
        if ns <= 0:
            continue
        acc = {c: 0.0 for c in chans}
        ass = 0.0
        for s in range(ns):
            t = s + dt
            ass += (-np.trace(Phi[t] @ tau[t, s] @ Phi[s] @ tau[s, t])).real
            ttPS = tau[s, s] - 0.5 * Iv
            # FS source equal-time leg: use -taugw (CONSISTENT with the -taugw bridging legs -- Fin sign fix),
            # tadpole-subtracted: ttFS = -taugw(s,s) + c_FS I,  c_FS = tr[Phi taugw]/tr[Phi]  (-> tr[Phi ttFS]=0)
            cFS = (np.trace(Phi[s] @ taugw[s, s]) / np.trace(Phi[s])).real
            ttFS = -taugw[s, s] + cFS * Iv
            src_out_PS = tau[s, t]
            src_out_FS = -taugw[s, t]
            sink_in_PS = tau[t, s]
            sink_in_FS = -taugw[t, s]
            for (X, Y) in chans:
                tt = ttFS if Y == "FS" else ttPS
                so = src_out_FS if Y == "FS" else src_out_PS
                si = sink_in_FS if X == "FS" else sink_in_PS
                acc[(X, Y)] += (-np.trace(Phi[s] @ tt @ Phi[s] @ so @ Phi[t] @ si)).real
        Css[dt] = ass / ns
        for c in chans:
            C[c][dt] = acc[c] / ns
    return C, Css, chans, twin


def main():
    C, Css, chans, twin = build()
    order = [("PS", "PS"), ("FS", "PS"), ("PS", "FS"), ("FS", "FS")]
    lab = {c: "%s<-%s" % c for c in order}
    # ABSOLUTE numerator effmass (Fin: read the state directly; physical constraint: sink sigma_PS/FS single
    # bilinear can reach 0.378 then 0.556 but NOT lighter than 0.378 and NOT the two-meson 0.756)
    eff = {}
    for c in order:
        with np.errstate(all="ignore"):
            eff[c] = np.log(np.abs(C[c][:-1] / C[c][1:]))
    print("# %s flavor triangles  <sigma_X(t) sigma_Y^2(0)>  ABSOLUTE effmass + R=C/Css (X<-Y)" % LTAG)
    print("# refs: m_sig=0.378/0.393; R plateaus to the relative single-meson OVERLAP (amplitude-vs-L test).")
    print("#  dt |  " + "   ".join("%-10s" % lab[c] for c in order) + " |  R(FS<-PS)   R(PS<-FS)   R(FS<-FS)")
    for dt in range(4, min(twin - 2, 22)):
        row = "   ".join("%+9.4f " % eff[c][dt] if np.isfinite(eff[c][dt]) else "   ---     " for c in order)
        rr = "   ".join("%+9.5f" % (C[c][dt] / Css[dt]) for c in [("FS", "PS"), ("PS", "FS"), ("FS", "FS")])
        print("#  %2d |  %s |  %s" % (dt, row, rr))


if __name__ == "__main__":
    main()
