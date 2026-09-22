#!/usr/bin/env python3
# sigma2_cross_trunc_check_claude.py -- reconcile with {1,1,1,1}: is the CROSS <sigma_00(t) sigma^2_00(0)> (C_12,
#   the actual <single-meson | sigma^2> overlap) zero and ROBUST to basis truncation, even while the sigma^2 AUTO
#   correlator C_22 ground collapses toward m_PS?  Loops NVKEEP over the free L1 config and prints, per NVKEEP:
#     R12 = |C_12| / |C_11|   (their normalized cross overlap; ~1e-7 = "exactly 0")   at a few dt,
#     C_22-ground vs C_11-ground effmass (to show the auto collapse is a SEPARATE effect).
#   If R12 stays ~0 at Nv=24,18,12,6 -> the exclusion <sigma^2|single-meson>=0 is truncation-robust (NOT a
#   completeness property); the C_22 shift is an auto-GEVP / degenerate-state artifact, not genuine m_PS overlap.
#   Run: ENS=free LREF=1 NVDIR=distill_Nv24 python3 sigma2_cross_trunc_check_claude.py

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

NVLIST = [int(x) for x in os.environ.get("NVLIST", "0,18,12,6").split(",")]  # 0 = complete (all modes)
KCFG = int(os.environ.get("KCFG", str(dc.KS[0])))
DTS = [int(x) for x in os.environ.get("DTS", "4,6,8,10").split(",")]


def eff(C):
    # log effmass of a real positive correlator vector C(dt)
    with np.errstate(all="ignore"):
        return np.log(C[:-1] / C[1:])


def main():
    dualf = dc.dual_areas_from_mesh().astype(float)
    wY = dualf * dc.Y00
    Pmap = MG.G.antipodal_map()
    print("# ENS=%s L=%d NVDIR=%s  k=%d   CROSS C_12=<sigma_00 sigma^2> vs AUTO C_11,C_22  (free)"
          % (dc.ENS.split("nu0")[0], dc.L, os.environ["NVDIR"], KCFG))
    print("# R12 = |C_12|/|C_11|  (~1e-7 => exclusion holds).  eff = -log ratio (single-cfg, no err).")
    hdr = "# Nv  | " + "  ".join("R12(dt%d)" % d for d in DTS) + " |  C11eff(dt8)  C22eff(dt8)"
    print(hdr)
    for nv in NVLIST:
        os.environ["NVKEEP"] = str(nv)                 # read fresh inside dc._read_peram_raw on each load
        C = MG.one_config(KCFG, dualf, wY, Pmap)       # (2,2,DTMAX), real parts
        C11 = np.abs(C[0, 0])
        C22 = np.abs(C[1, 1])
        C12 = np.abs(C[0, 1])
        r12 = [C12[d] / C11[d] if C11[d] != 0 else np.nan for d in DTS]
        e11 = eff(C11)
        e22 = eff(C22)
        lab = "all" if nv == 0 else "%d" % nv
        print("# %4s | %s | %10.4f  %10.4f"
              % (lab, "  ".join("%9.2e" % x for x in r12), e11[8], e22[8]))


if __name__ == "__main__":
    main()
