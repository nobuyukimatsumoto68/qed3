#!/usr/bin/env python3
# contact_test_claude.py
#   >>> SUPERSEDED / ARTIFACT (2026-09-17).  This script's "St"/"S+St" modes call AblkSt (the FS Stilde
#   >>> furnished leg = -taugw), which is the tau_gw ARTIFACT and is now DISABLED at the call sites
#   >>> (fs_gevp_point FS_STILDE_FURNISH_ENABLED=False).  The FS Stilde leg collapses to the S-part (tau) by
#   >>> GW -- see fs_furnishing_derivation_claude.md.  So the S/Stilde split explored here is refuted; the
#   >>> "St" number is not physical.  Kept only as the investigation record; use the "S" (AblkS=tau) mode.
#   Does the one-meson (m_sigma=0.378) in the sigma^2_00 correlator come from the equal-time CONTACT term?
#   Build the full connected sigma^2_00 (op 0) diagonal correlator two ways and compare effmass:
#     (a) contact OUT  -- AblkS as shipped (equal-time block = tau_eq - 1/2 I = tilde_tau)  [current code]
#     (b) contact IN   -- add the 1/2 I back at equal time
#   If (a) removes the 0.378 and (b) shows it, the one-meson is a pure contact artifact (NM's claim).
#   Free limit, single exact config.  Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 contact_test_claude.py

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
import fs_gevp_point_claude as G

DTMAX = int(os.environ.get("DTMAX", "24"))


def sigma2_diag(k, contact_in, legs):
    # legs = "S", "St", or "S+St"
    AblkS0, AblkSt, twin, nsite, U, tau = G.make_config(k)
    NS = G.NS
    # identity block on (site,spin) to add the contact back
    Iblk = np.zeros((nsite, NS, nsite, NS))
    for s in range(nsite):
        for a in range(NS):
            Iblk[s, a, s, a] = 1.0

    def AblkS(ta, tb):
        A = AblkS0(ta, tb)
        if contact_in and ta == tb:
            A = A + 0.5 * Iblk                       # undo the -1/2 contact removal
        return A

    which = {"S": [AblkS], "St": [AblkSt], "S+St": [AblkS, AblkSt]}[legs]
    dual = dc.dual_areas_from_mesh().astype(float)
    wY = dual * dc.Y00
    Pmap = G.antipodal_map()
    C = np.full(DTMAX, np.nan)
    for dt in range(DTMAX):
        s0s = [s for s in range(twin) if s + dt < twin]
        if not s0s:
            continue
        s0s = np.array(s0s)
        vt = [dt, dt, 0, 0]
        offs = set((vt[va], vt[vb]) for va in range(4) for vb in range(4))
        vspec = G.op_vspec(0, ('i', 'j'), dual, wY) + G.op_vspec(0, ('k', 'l'), dual, wY)
        v = 0.0
        for blk in which:
            bA = {o: np.array([blk(s + o[0], s + o[1]) for s in s0s]) for o in offs}
            for cyc in G.PERMS:
                v += G.perm_contrib_folded(cyc, vt, bA, vspec, Pmap)
        C[dt] = (v.real) / len(s0s)
    return C


def eff(C):
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


def main():
    k = dc.KS[0]
    cS = sigma2_diag(k, contact_in=False, legs="S")            # S leg only, contact out (PS-like)
    cSt = sigma2_diag(k, contact_in=False, legs="St")          # Stilde leg only, contact out
    cFS = sigma2_diag(k, contact_in=False, legs="S+St")        # FS (S+Stilde), contact out (current 3x3)
    mS = eff(cS)
    mSt = eff(cSt)
    mFS = eff(cFS)
    print("# sigma^2_00 diagonal effmass by furnishing leg (all contact OUT)")
    print("# refs: one-meson m_sigma=0.378 ; (2,2)=0.556 ; two-meson=0.756")
    print("#  dt |  m_S(PS)   m_St     m_FS(current)")
    for dt in range(1, DTMAX - 1):
        r = []
        for m in (mS, mSt, mFS):
            r.append("%7.4f" % m[dt] if np.isfinite(m[dt]) else "  ---  ")
        print("#  %2d | %s" % (dt, "  ".join(r)))


if __name__ == "__main__":
    main()
