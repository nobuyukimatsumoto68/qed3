#!/usr/bin/env python3
# hankel_rebase_runscan_claude.py
#   Controlled one-parameter-at-a-time scans of the block-Hankel + rebase knobs, using the
#   cached store from hankel_rebase_scan_claude.py.  Prints an effmass table + a summary
#   (window-averaged m0,m1 with jackknife errors and the m1 slope) for each configuration,
#   grouped so the TREND of each knob is visible.  Does not modify any existing driver.

import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("ENS", "Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L1_hb1.000000")
os.environ.setdefault("NVDIR", "distill_Nv24")
import sys
sys.path.insert(0, ".")
import numpy as np
import hankel_rebase_scan_claude as hs

store, twin, ncfg, tag = hs.build_store()


def report(name, nooa, offsets, stages, t0, w0, w1, nkeep_show=2):
    em_c, em_err, ems, nbin = hs.eval_config(store, twin, nooa, offsets, stages, t0)
    tmax = em_c.shape[0]
    nkeep = em_c.shape[1]
    print("\n=== %s ===" % name)
    print("    nooa=%d offsets=%s stages=%s T0=%d tmax=%d nkeep=%d" %
          (nooa, offsets, stages, t0, tmax, nkeep))
    hdr = "  t  | " + "  ".join("m%d(err)" % n for n in range(min(nkeep, nkeep_show)))
    print(hdr)
    for t in range(t0 + 1, min(tmax, 22)):
        cells = []
        for n in range(min(nkeep, nkeep_show)):
            v = em_c[t, n]
            e = em_err[t, n]
            cells.append("%6.4f(%.4f)" % (v, e) if np.isfinite(v) else "   ---    ")
        print("  %2d | %s" % (t, "  ".join(cells)))
    # summary: window mass + slope
    m0, e0 = hs.window_mass(em_c, ems, nbin, 0, w0[0], w0[1])
    s0, se0 = hs.window_slope(em_c, ems, nbin, 0, w0[0], w0[1])
    line = "  SUMMARY m0=%.4f(%.4f) [win %d-%d slope %+.4f(%.4f)]" % (m0, e0, w0[0], w0[1], s0, se0)
    if nkeep >= 2:
        m1, e1 = hs.window_mass(em_c, ems, nbin, 1, w1[0], w1[1])
        s1, se1 = hs.window_slope(em_c, ems, nbin, 1, w1[0], w1[1])
        tilt = "TILT!" if abs(s1) > 2 * se1 and abs(s1) > 0.008 else "flat"
        line += "\n           m1=%.4f(%.4f) [win %d-%d slope %+.4f(%.4f) -> %s]" % (
            m1, e1, w1[0], w1[1], s1, se1, tilt)
    print(line)
    return dict(name=name, nooa=nooa, offsets=offsets, stages=stages, t0=t0, em_c=em_c,
                em_err=em_err, tmax=tmax,
                m0=m0, e0=e0, s0=s0, se0=se0,
                m1=(m1 if nkeep >= 2 else np.nan), e1=(e1 if nkeep >= 2 else np.nan),
                s1=(s1 if nkeep >= 2 else np.nan), se1=(se1 if nkeep >= 2 else np.nan))


results = []


def section(title):
    print("\n" + "#" * 78)
    print("# " + title)
    print("#" * 78)


# ----------------------------------------------------------------------------
# A) BASELINE reproduce + T0 trend  (Dt=1,2 = offsets 0,1,2 ; 2-op ; rebase@5)
# ----------------------------------------------------------------------------
section("A) T0 trend  (Dt=1,2 [0,1,2], 2-op, single rebase 2@t=5)")
for t0 in [2, 3, 4]:
    r = report("A_T0=%d_Dt12_2op" % t0, 1, [0, 1, 2], [(5, 2)], t0, (5, 12), (6, 12))
    results.append(r)

# ----------------------------------------------------------------------------
# B) Hankel-shape trend (fix T0=3, 2-op, rebase 2@t=5)
# ----------------------------------------------------------------------------
section("B) Hankel-shape trend  (2-op, T0=3, single rebase 2@t=5)")
shapes = [
    ("B_Dt1_[0,1]", [0, 1], (6, 12)),
    ("B_Dt12_[0,1,2]", [0, 1, 2], (6, 12)),
    ("B_Dt24_[0,2,4]", [0, 2, 4], (6, 12)),
    ("B_Dt3_[0,3]", [0, 3], (6, 12)),
    ("B_Dt123_[0,1,2,3]", [0, 1, 2, 3], (6, 12)),
    ("B_Dt12345_[0..5]", [0, 1, 2, 3, 4, 5], (6, 12)),
    ("B_Dt48_[0,4,8]_REJECTcand", [0, 4, 8], (6, 12)),
]
for nm, off, w1 in shapes:
    r = report(nm, 1, off, [(5, 2)], 3, (5, 12), w1)
    results.append(r)

# ----------------------------------------------------------------------------
# C) 2-op vs 3-op base trend (fix Dt=1,2, T0=3, rebase@5)
# ----------------------------------------------------------------------------
section("C) 2-op vs 3-op base  (Dt=1,2 [0,1,2], T0=3, rebase@5)")
r = report("C_2op_nkeep2", 1, [0, 1, 2], [(5, 2)], 3, (5, 12), (6, 12), nkeep_show=2)
results.append(r)
r = report("C_3op_nkeep2", 0, [0, 1, 2], [(5, 2)], 3, (5, 12), (6, 12), nkeep_show=2)
results.append(r)
r = report("C_3op_nkeep3", 0, [0, 1, 2], [(5, 3)], 3, (5, 12), (6, 12), nkeep_show=3)
results.append(r)

# ----------------------------------------------------------------------------
# D) rebase-time trend for 3-op (Dt=1,2 [0,1,2], T0=3, nkeep2)
# ----------------------------------------------------------------------------
section("D) rebase-time trend (3-op nkeep2, Dt=1,2 [0,1,2], T0=3)")
for rt in [4, 5, 7, 9]:
    r = report("D_3op_reb%d" % rt, 0, [0, 1, 2], [(rt, 2)], 3, (5, 12), (6, 12))
    results.append(r)

# ----------------------------------------------------------------------------
# E) staged (multiple) rebase (3-op, fine ladder)
# ----------------------------------------------------------------------------
section("E) staged rebase  (3-op, fine ladder Dt=1,2,3 [0,1,2,3], T0=3)")
r = report("E_3op_single_2@5", 0, [0, 1, 2, 3], [(5, 2)], 3, (5, 12), (6, 12))
results.append(r)
r = report("E_3op_staged_4@4_2@8", 0, [0, 1, 2, 3], [(4, 4), (8, 2)], 3, (5, 12), (6, 12))
results.append(r)
r = report("E_3op_staged_6@4_2@9", 0, [0, 1, 2, 3, 4, 5], [(4, 6), (9, 2)], 3, (5, 12), (6, 12))
results.append(r)

# ----------------------------------------------------------------------------
# ranked summary
# ----------------------------------------------------------------------------
print("\n" + "#" * 78)
print("# RANKED SUMMARY (all configs)")
print("#" * 78)
print("%-30s %10s %10s %12s %6s" % ("name", "m0", "m1", "m1-slope", "flat?"))
for r in results:
    flat = "flat" if (np.isfinite(r["s1"]) and abs(r["s1"]) <= 2 * r["se1"]) or \
        (np.isfinite(r["s1"]) and abs(r["s1"]) <= 0.008) else "TILT"
    print("%-30s %5.4f(%.4f) %5.4f(%.4f) %+7.4f(%.4f) %6s" %
          (r["name"], r["m0"], r["e0"], r["m1"], r["e1"], r["s1"], r["se1"], flat))

import pickle
with open(hs.SCRATCH + "/scan_results_claude.pkl", "wb") as f:
    pickle.dump(results, f)
print("\n# results pickled -> %s/scan_results_claude.pkl" % hs.SCRATCH)
