#!/usr/bin/env python3
# glue_f2_v2_compare_claude.py -- Chunks 4 and 6 of glue_f2_v2_shapes_impl_plan_claude.md.
#
# Two modes, both driven by run_glue_gevp_f2_v2_claude.sh:
#
#   gate <prod.dat> <gate.dat>
#       VALIDATION GATE. The p=2 sub-block of the power-extended basis IS the production F^2 basis,
#       so the two GEVP outputs must agree. The driver computes obs for shapes 0..6 with op_pow(p=2)
#       = g=th; g*=th, i.e. the same single multiply and the same accumulation order as the old
#       op(squared=true), and the correlator entries for a,b <= 6 are the same products summed over
#       the same t -- so BIT-IDENTICAL is the expectation, not merely "close".
#       exit 0 = pass, 1 = fail (run_glue_gevp_f2_v2_claude.sh skips PHASE B on a fail).
#
#   verdict <tag> "<nops2 list>" "<vacsub list>"
#       Fit every state of every full-basis variant over the F^2 window and tabulate against the
#       production 0++, so the state carrying the 0++ can be identified (the production bookkeeping
#       nops2=2 / state 0 does NOT survive the enlarged l=0 block) and the error ratio read off.
#
# The constant fit is delegated to fit_perm_claude.py (diagonal-chi2 center, correlation-aware
# jackknife error), the same estimator the production plots use.
import subprocess
import sys
import numpy as np

# F^2 0++ constant-fit window, matching FITS["f2"] in plot_glue_gevp_claude.py
TLO = 0.2
THI = 0.4

# gate tolerances on the max relative difference between the two .dat tables
TOL_EXACT = 1.0e-12   # at or below: bit-identical, as expected
TOL_FAIL = 1.0e-06    # above: fail


def fit_state(jk, state):
    # returns (M, err) for one state, via fit_perm_claude.py
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(TLO), str(THI), str(state)],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if ("state %d:" % state) in line:
            seg = line.split("M=")[-1]
            tok = seg.strip().split()[0]
            val = float(tok.split("(")[0])
            err = float(tok.split("(")[1].rstrip(")"))
            return val, err
    return None, None


def nstates_of(jk):
    with open(jk) as f:
        h = f.readline().split()
    return int(h[3])


def mode_gate(fn_prod, fn_gate):
    a = np.loadtxt(fn_prod)
    b = np.loadtxt(fn_gate)
    print("# GATE: %s  vs  %s" % (fn_prod, fn_gate))
    if a.shape != b.shape:
        print("  FAIL: shape mismatch %s vs %s" % (a.shape, b.shape))
        return 1
    # NaN-aware: the cosh effective mass is legitimately NaN once the correlator turns over (here
    # rows t >= 2.0), and nan-nan would propagate into every comparison. Require the SAME NaN
    # pattern, then compare only the finite entries.
    na = np.isnan(a)
    nb = np.isnan(b)
    if not np.array_equal(na, nb):
        print("  FAIL: NaN pattern differs (prod %d, gate %d NaNs)" % (na.sum(), nb.sum()))
        return 1
    fin = ~na
    print("  table shape       : %s  (col 0 = t, cols 3,4 = state 0 = the 0++)" % (a.shape,))
    print("  NaN entries       : %d, identical pattern in both" % na.sum())
    den = np.maximum(np.abs(a[fin]), 1.0e-300)
    rel = np.abs(a[fin] - b[fin]) / den
    print("  max |rel diff|    : %.3e   (over %d finite entries)" % (rel.max(), fin.sum()))
    print("  max |abs diff|    : %.3e" % np.abs(a[fin] - b[fin]).max())
    print("  0++ (state 0) t=%.2f : prod %.10e  gate %.10e" % (a[0, 0], a[0, 3], b[0, 3]))
    if rel.max() <= TOL_EXACT:
        print("  PASS (bit-identical to %.0e)" % TOL_EXACT)
        return 0
    if rel.max() <= TOL_FAIL:
        print("  PASS WITH WARNING: agreement is only %.3e, not bit-identical." % rel.max())
        print("  Expected 0 here -- investigate before trusting PHASE B.")
        return 0
    print("  FAIL: %.3e exceeds %.0e" % (rel.max(), TOL_FAIL))
    return 1


def mode_verdict(tag, nops2_list, vacsub_list):
    ref_jk = "gevp_f2v2_prod_%s_jk_claude.dat" % tag
    m_ref, e_ref = fit_state(ref_jk, 0)
    print("")
    print("=== F^2 0++ verdict, %s, constant fit over t in [%.1f, %.1f] ===" % (tag, TLO, THI))
    if m_ref is None:
        print("could not fit the production reference %s -- aborting table" % ref_jk)
        return
    print("production (7 ops, p=2, nops2=2, vacsub=0): state 0  M = %.4f(%.4f)" % (m_ref, e_ref))
    print("")
    print("%-18s %-7s %-7s %-9s %-9s %-8s" % ("variant", "nops2", "vacsub", "state", "M(err)", "err/ref"))
    print("-" * 68)
    for vs in vacsub_list:
        for n2 in nops2_list:
            name = "full_n%s_v%s" % (n2, vs)
            jk = "gevp_f2v2_%s_%s_jk_claude.dat" % (name, tag)
            try:
                ns = nstates_of(jk)
            except OSError:
                print("%-18s %-7s %-7s  (missing: %s)" % (name, n2, vs, jk))
                continue
            best = None
            for s in range(ns):
                m, e = fit_state(jk, s)
                if m is None:
                    continue
                ratio = e / e_ref if e_ref > 0 else float("nan")
                mark = ""
                if best is None or abs(m - m_ref) < abs(best[1] - m_ref):
                    best = (s, m)
                print("%-18s %-7s %-7s %-9d %.4f(%.4f) %-8.2f"
                      % (name, n2, vs, s, m, e, ratio))
            if best is not None:
                print("%-18s %-7s %-7s  -> state %d is closest to the production 0++"
                      % ("", "", "", best[0]))
            print("-" * 68)
    print("")
    print("Decision rule (plan section 6): keep the enlargement only if the 0++ error drops by more")
    print("than its own jackknife uncertainty. err/ref < 1 is necessary, not sufficient.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: glue_f2_v2_compare_claude.py gate <prod.dat> <gate.dat>")
        print("       glue_f2_v2_compare_claude.py verdict <tag> \"<nops2 list>\" \"<vacsub list>\"")
        sys.exit(2)
    if sys.argv[1] == "gate":
        sys.exit(mode_gate(sys.argv[2], sys.argv[3]))
    if sys.argv[1] == "verdict":
        mode_verdict(sys.argv[2], sys.argv[3].split(), sys.argv[4].split())
        sys.exit(0)
    print("unknown mode: %s" % sys.argv[1])
    sys.exit(2)
