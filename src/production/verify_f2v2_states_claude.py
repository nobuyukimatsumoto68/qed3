# verify_f2v2_states_claude.py -- confirm the physical 0++ is STATE 2 in every gevp_f2 jk dump
# (F^2-v2 / F^4 basis, nops2=4, vacsub=0).  Fits each of states 0..3 over the 0++ window [0.2,0.4]
# and checks that state 2 is in a physical range and is the intended 0++ (not the ~0 vacuum, not nan,
# not the higher ~6 state).  Flags any ensemble where the ordering breaks so the plot column choice
# (cols 7,8 = state 2) is not silently wrong.
import glob
import re
import fit_perm_claude as fp
import numpy as np

TLO, THI = 0.2, 0.4
PHYS_LO, PHYS_HI = 1.0, 5.5   # a physical 0++ a_t m sits here; vacuum ~0, higher state ~6

def fit_state(jk, s):
    a, at = fp.load(jk)
    ti = [i for i in range(a.shape[1]) if TLO - 1e-9 <= (i + 1) * at <= THI + 1e-9]
    Mb, M, sig, chi2, dof = fp.fit_state(a, ti, s)
    return M, sig

rows = []
for jk in sorted(glob.glob("gevp_f2_Nf*_L*_jk_claude.dat")):
    tag = re.sub(r"^gevp_f2_|_jk_claude\.dat$", "", jk)
    ms = []
    for s in range(4):
        try:
            M, sig = fit_state(jk, s)
        except Exception:
            M, sig = float("nan"), float("nan")
        ms.append(M)
    s2 = ms[2]
    ok = (not np.isnan(s2)) and PHYS_LO <= s2 <= PHYS_HI
    # is state 2 actually the one closest to physical, or did a different state take that role?
    phys_states = [s for s in range(4) if not np.isnan(ms[s]) and PHYS_LO <= ms[s] <= PHYS_HI]
    rows.append((tag, ms, ok, phys_states))

print("%-22s %9s %9s %9s %9s  verdict" % ("ensemble", "state0", "state1", "state2", "state3"))
print("-" * 78)
nbad = 0
for tag, ms, ok, phys in rows:
    cells = " ".join("%9.3f" % m if not np.isnan(m) else "      nan" for m in ms)
    if ok and phys == [2]:
        v = "OK (0++=s2)"
    elif ok:
        v = "s2 ok but ALSO physical in states %s" % phys
    else:
        v = "*** BAD: s2=%.3f not physical; physical in %s ***" % (ms[2], phys)
        nbad += 1
    print("%-22s %s  %s" % (tag, cells, v))
print("-" * 78)
print("%d ensembles, %d with a state-2 problem" % (len(rows), nbad))
