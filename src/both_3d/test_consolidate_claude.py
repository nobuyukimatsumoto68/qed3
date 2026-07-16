#!/usr/bin/env python3
# Validate operator consolidation: sum the orbits within each shape (equal weight per orbit) -> 5
# shape-operators, and check the per-m GEVP ground matches the full n_orbits basis. L2 F l=1, Nf2 g8.
# Replicates the analysis fold: matCA[k] = sym(Fc[k+1]) + sym(Fc[k+2]); GEVP(matCA[ti],matCA[ti+1]).
import glob
import re
import numpy as np
import h5py

DD = "data_Nf2_gsq8.000000at0.200000nu01.000000nt128L2"
AT = 0.2
BINSIZE = 50
NOB = 14
NLM = 8
RTOL = 1e-8
TREBASE = 1
# L2 orbit -> shape sizes: triangle2, rect2, fig8 5, three-tri3, star2
SHAPE_SIZES = [2, 2, 5, 3, 2]
SHNAME = ["triangle", "rect", "fig8", "three-tri", "star"]

# consolidation matrix T (n_shapes x n_orbits), equal weight (1) per orbit
T = np.zeros((5, NOB))
o = 0
for s, n in enumerate(SHAPE_SIZES):
    for _ in range(n):
        T[s, o] = 1.0
        o += 1

fs = glob.glob(f"{DD}/glue_msm_shapes.[0-9]*.h5")
ks = sorted(int(re.search(r'\.(\d+)\.h5$', f).group(1)) for f in fs)
Nc = len(ks)
nbins = Nc // BINSIZE
N = nbins * BINSIZE
nsep = 17
Nthalf = nsep - 2

def inv_sqrt(M):
    w, V = np.linalg.eigh(M)
    wmax = w.max()
    idv = np.array([1.0/np.sqrt(x) if (x > RTOL*wmax and x > 0) else 0.0 for x in w])
    return V @ np.diag(idv) @ V.T

def sym(X):
    return 0.5*(X + X.T)

def ground_effmass(C_by_sep, ti):
    # C_by_sep[sep] = n x n correlator at separation sep (sep=0..16); returns ground effmass at output t=ti*at
    matCA = [sym(C_by_sep[k+1]) + sym(C_by_sep[k+2]) for k in range(Nthalf)]  # k=0..14
    si0 = inv_sqrt(matCA[TREBASE])
    M0 = si0 @ matCA[TREBASE+1] @ si0
    _, Vre = np.linalg.eigh(sym(M0))
    si = inv_sqrt(matCA[ti])
    M = si @ matCA[ti+1] @ si
    Mrot = Vre.T @ sym(M) @ Vre
    ev = np.linalg.eigvalsh(Mrot)
    lam = ev[-1]   # ground = largest
    return -np.log(lam)/AT if lam > 0 else np.nan

# read the l=1 orbit correlators (ilm=0,1,2 = m=-1,0,1), full-sample average
# blk[dt][ilm*NOB*NOB + a*NOB + b] = C(orbit a, orbit b) channel ilm
acc = np.zeros((3, nsep, NOB, NOB))
for k in ks[:N]:
    with h5py.File(f"{DD}/glue_msm_shapes.{k}.h5", "r") as h:
        blk = h["F_corr_blk"][:]
    for ilm in range(3):
        acc[ilm] += blk[:, ilm*NOB*NOB:(ilm+1)*NOB*NOB].reshape(nsep, NOB, NOB)
acc /= N

print("Per-m GROUND effmass: full 14-orbit vs consolidated 5-shape (Nf2 g8 L2 F l=1, full sample)")
print(f"{'t':>5} | {'m=-1 full':>10} {'cons':>7} | {'m=0 full':>10} {'cons':>7} | {'m=+1 full':>10} {'cons':>7}")
for ti in [1, 2, 3, 4]:   # t = 0.2..0.8
    row = f"{ti*AT:5.1f} |"
    for ilm in range(3):
        Cfull = [acc[ilm][sep] for sep in range(nsep)]
        Ccons = [T @ acc[ilm][sep] @ T.T for sep in range(nsep)]
        gf = ground_effmass(Cfull, ti)
        gc = ground_effmass(Ccons, ti)
        row += f" {gf:10.4f} {gc:7.4f} |"
    print(row)
