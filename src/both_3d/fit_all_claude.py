#!/usr/bin/env python3
# Correlated constant-fit driver: loop the reliable ensembles, fit the plateau for F (l=1, all 3
# degenerate levels + first NFIT_F t-points, simultaneously) and F^2 (l=0, state 0, first NFIT_F2),
# using the jackknife covariance (fit_plateau_claude.fit_const). Emit a summary table and per-Nf
# trend data files (M vs gsq) for gnuplot.
#
# Fit ranges are retunable here (NFIT_F / NFIT_F2 / L2 overrides). Singular ensembles (too few bins
# for the covariance, i.e. the thin ~320-cfg couplings) are reported but excluded from the trend.
import numpy as np
from fit_plateau_claude import load_jk, fit_const

NFIT_F  = 8     # F l=1: first 8 t-points x 3 levels = 24 correlated points
NFIT_F2 = 3     # F^2 0++: signal only at small t -> first 3 points (retune)
NFIT_F_L2 = 4   # L2 F: short plateau -> first 4 points (retune)

# reliable list (gsq 2.2, 2.5 excluded). (Nf, gsq, L)
ENS = [
    (2,"8.000000",1),(4,"8.000000",1),(6,"8.000000",1),
    (2,"8.000000",2),(4,"8.000000",2),(6,"8.000000",2),
    (2,"1.000000",1),(2,"2.000000",1),(2,"2.400000",1),(2,"4.000000",1),(2,"12.000000",1),
    (4,"1.000000",1),(4,"2.000000",1),(4,"4.000000",1),(4,"12.000000",1),
    (6,"1.000000",1),(6,"2.000000",1),(6,"2.400000",1),(6,"4.000000",1),(6,"12.000000",1),
]

def do(tag, kind, states, nfit):
    fn = f"jk{kind}_{tag}_claude.dat"
    try:
        jk, at = load_jk(fn)
    except Exception:
        return None
    ntpts = jk.shape[1]
    ti = list(range(min(nfit, ntpts)))
    return fit_const(jk, ti, states)

rows = []
trendF = {}     # (Nf,L) -> list of (gsq, M, err, chi2dof)
trendF2 = {}
print(f"{'ensemble':22s} | {'F l=1 (24pt)':>26s} | {'F^2 0++':>24s}")
for nf,g,L in ENS:
    tag = f"Nf{nf}_g{g}_L{L}"
    nfitF = NFIT_F_L2 if L>=2 else NFIT_F
    rF  = do(tag, "F",  [0,1,2], nfitF)
    rF2 = do(tag, "F2", [0],     NFIT_F2)
    def s(r):
        if r is None: return "   --   "
        if r["singular"]: return f"singular(nb={r['nbins']})"
        return f"{r['M']:.3f}({r['err']:.3f}) X2/d={r['chi2']/r['dof']:.1f}"
    print(f"Nf{nf} g{g[:5]:5s} L{L}      | {s(rF):>26s} | {s(rF2):>24s}")
    gq = float(g)
    if rF and not rF["singular"]:
        trendF.setdefault((nf,L),[]).append((gq,rF["M"],rF["err"],rF["chi2"]/rF["dof"]))
    if rF2 and not rF2["singular"]:
        trendF2.setdefault((nf,L),[]).append((gq,rF2["M"],rF2["err"],rF2["chi2"]/rF2["dof"]))

# write trend data files: gsq M err chi2dof, sorted by gsq
for kind,tr in [("F",trendF),("F2",trendF2)]:
    for (nf,L),lst in tr.items():
        lst.sort()
        out = f"trend_{kind}_Nf{nf}_L{L}_claude.dat"
        with open(out,"w") as f:
            f.write("# gsq  M  err  chi2/dof\n")
            for gq,M,e,c in lst:
                f.write(f"{gq:.4f} {M:.6f} {e:.6f} {c:.3f}\n")
        print("wrote", out)
