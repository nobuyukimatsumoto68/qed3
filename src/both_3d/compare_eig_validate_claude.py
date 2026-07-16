#!/usr/bin/env python3
# compare_eig_validate_claude.py
# _claude: validate the IRL (eig_lanczos) against the dense geev reference (eig_wmass_val) by diffing the
# lowest Nk eigenvalues of A = (D_ov+m)^dag (D_ov+m).  Both .dat files list "i  eval  sqrt(eval)", ascending.
# PASS if every one of the lowest Nk IRL eigenvalues matches the dense reference to the tolerances below.
# Usage: python3 compare_eig_validate_claude.py <irl.dat> <dense.dat> [Nk]

import sys
import numpy as np


def load(path):
    vals = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            vals.append(float(parts[1]))
    return np.array(sorted(vals))


irl_path = sys.argv[1]
dense_path = sys.argv[2]
Nk = int(sys.argv[3]) if len(sys.argv) > 3 else 20

irl = load(irl_path)
dense = load(dense_path)

n = min(Nk, len(irl), len(dense))
print(f"# comparing lowest {n} eigenvalues  (IRL has {len(irl)}, dense has {len(dense)})")
print(f"# {'i':>3} {'IRL':>22} {'dense':>22} {'abs_diff':>12} {'rel_diff':>12}")

ATOL = 1.0e-6
RTOL = 1.0e-5
max_abs = 0.0
max_rel = 0.0
ok = True
for i in range(n):
    a = irl[i]
    b = dense[i]
    ad = abs(a - b)
    rd = ad / max(abs(b), 1.0e-30)
    max_abs = max(max_abs, ad)
    max_rel = max(max_rel, rd)
    flag = "" if (ad < ATOL or rd < RTOL) else "  <-- MISMATCH"
    if flag:
        ok = False
    print(f"  {i:>3} {a:>22.14e} {b:>22.14e} {ad:>12.3e} {rd:>12.3e}{flag}")

print(f"# max abs_diff = {max_abs:.3e}   max rel_diff = {max_rel:.3e}")
print(f"# tolerances: atol={ATOL:.1e} rtol={RTOL:.1e} (pass if EITHER holds per mode)")
if ok:
    print("# RESULT: PASS -- IRL reproduces the dense low spectrum")
    sys.exit(0)
else:
    print("# RESULT: FAIL -- see MISMATCH rows above")
    sys.exit(1)
