#!/usr/bin/env python3
# one_point_claude.py
# Calculate the equal-time ONE-POINT functions of the scalar operator from the data, to see how far
# they sit beyond the analytic GW contact (1/2).  Building blocks (PS legs = tau), per window-time s:
#   D_S(s)   = Tr[ Phi(s) ( tau(s,s) - c I ) ]                     (single-sigma tadpole / condensate)
#   D'_S(s)  = Tr[ Phi(s)(tau(s,s)-cI) Phi(s)(tau(s,s)-cI) ]       (the D' one-point)
#   O(s)     = 2 ( D_S(s)^2 + D'_S(s) )                            (composite sigma^2 one-point)
# c = 0 (raw) and c = 1/2 (contact-subtracted).  Report ensemble means (config jackknife) + the
# purely-analytic contact contribution c*Tr[Phi] for reference.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys
sys.path.insert(0, ".")
import numpy as np
import distill_contract_claude as dc


def onepoints(c):
    dual = dc.dual_areas_from_mesh()
    w00 = np.repeat(dual, dc.NS) * dc.Y00
    DS = []            # per-config window-mean D_S
    DpS = []           # per-config window-mean D'_S
    O = []             # per-config window-mean O = 2(D_S^2 + D'_S)
    trPhi = []         # per-config window-mean Tr[Phi] (for the analytic contact c*Tr[Phi])
    for k in dc.KS:
        V, tau, taugw, tsrc0, twin = dc.load_peram(k)
        Iv = np.eye(tau.shape[-1])
        ds = 0.0
        dps = 0.0
        oo = 0.0
        tp = 0.0
        for s in range(twin):
            Vt = V[tsrc0 + s].T
            Phi = Vt.conj().T @ (w00[:, None] * Vt)
            tss = tau[s, s] - c * Iv
            d = np.trace(Phi @ tss)
            dp = np.trace(Phi @ tss @ Phi @ tss)
            ds += d.real
            dps += dp.real
            oo += (2.0 * (d ** 2 + dp)).real
            tp += np.trace(Phi).real
        DS.append(ds / twin)
        DpS.append(dps / twin)
        O.append(oo / twin)
        trPhi.append(tp / twin)
    return np.array(DS), np.array(DpS), np.array(O), np.array(trPhi)


def jk(x):
    n = x.shape[0]
    samp = np.array([np.delete(x, i, 0).mean(0) for i in range(n)])
    return samp.mean(0), np.sqrt((n - 1) * np.mean((samp - samp.mean(0)) ** 2, 0))


def main():
    tag = dc.ENS.split("nu0")[0]
    print("# ENS=%s  ncfg=%d  scalar one-point functions from data" % (tag, len(dc.KS)))
    for c in (0.0, 0.5):
        DS, DpS, O, trPhi = onepoints(c)
        mDS, eDS = jk(DS)
        mDp, eDp = jk(DpS)
        mO, eO = jk(O)
        contact = c * trPhi.mean()          # purely-analytic contact contribution to D_S
        print("\n=== contact subtraction c = %.2f ===" % c)
        print("  <D_S>   = %+.6f (+-%.6f)      [analytic contact c*Tr[Phi] = %+.6f]" % (mDS, eDS, contact))
        print("  <D'_S>  = %+.6f (+-%.6f)" % (mDp, eDp))
        print("  <O>=<2(D_S^2+D'_S)> = %+.6f (+-%.6f)" % (mO, eO))
        print("  std(D_S) over cfg   = %.3e   (fluctuation beyond the mean)" % DS.std(ddof=1))
        print("  std(D'_S) over cfg  = %.3e" % DpS.std(ddof=1))


if __name__ == "__main__":
    main()
