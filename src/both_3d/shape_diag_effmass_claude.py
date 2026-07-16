#!/usr/bin/env python3
# Per-shape, per-m DIAGONAL effective mass for the LINEAR F (l=1), to scrutinize the rank-1 claim.
# Uses the stored per-(l,m) blocks directly (no GEVP): for shape a and channel ilm, the diagonal
# self-correlator D = F_corr_blk[dt][ilm*nob^2 + a*nob + a]. Fold matches the analysis (windowed):
#   folded[k] = D[k+1] + D[k+2]   (k=0..Nthalf-1);  effmass(ti) = -log(folded[ti+1]/folded[ti])/at,
#   output t = ti*at,  ti = 1..tcut-1.   Jackknife over bins (binsize=50) for the error.
# If linear F is rank-1, all shapes AND all m collapse onto one curve.
import glob
import re
import numpy as np
import h5py

DD = "data_Nf4_gsq2.000000at0.200000nu01.000000nt128L1"
AT = 0.2
BINSIZE = 50
TCUT = 14
NOB = 5
NLM = 8
SHAPES = {0:"triangle", 1:"rect", 2:"fig8"}
MCH = {0:"m=-1", 1:"m=0", 2:"m=+1"}   # ilm 0,1,2 -> l=1 m=-1,0,1

fs = glob.glob(f"{DD}/glue_msm_shapes.[0-9]*.h5")
ks = sorted(int(re.search(r'\.(\d+)\.h5$', f).group(1)) for f in fs)
Nc = len(ks)
nbins = Nc // BINSIZE
N = nbins * BINSIZE
print(f"# Nc={Nc} nbins={nbins} N={N}")

# read all diagonal correlators D[config][a][ilm][sep]  (a in 0..2, ilm in 0..2, sep 0..16)
nsep = 17
D = np.zeros((N, 3, 3, nsep))
for ci in range(N):
    with h5py.File(f"{DD}/glue_msm_shapes.{ks[ci]}.h5", "r") as h:
        blk = h["F_corr_blk"][:]                 # (nsep, NLM*NOB*NOB)
    for a in range(3):
        for ilm in range(3):
            D[ci, a, ilm, :] = blk[:, ilm*NOB*NOB + a*NOB + a]
    if ci % 1000 == 0:
        print(f"# read {ci}/{N}")

# per-bin sums for jackknife
Nthalf = nsep - 2   # 15
tis = list(range(1, TCUT))          # ti = 1..13 -> t = 0.2..2.6

def effmass_from_D(Davg):
    # Davg[a,ilm,sep]; folded[k]=D[k+1]+D[k+2]; eff(ti)=-log(folded[ti+1]/folded[ti])/at
    out = np.full((3,3,len(tis)), np.nan)
    for a in range(3):
        for ilm in range(3):
            d = Davg[a,ilm]
            folded = np.array([d[k+1]+d[k+2] for k in range(Nthalf)])
            for n,ti in enumerate(tis):
                if folded[ti] > 0 and folded[ti+1] > 0:
                    out[a,ilm,n] = -np.log(folded[ti+1]/folded[ti]) / AT
    return out

# full-sample and jackknife (leave-one-bin-out)
binsum = D[:N].reshape(nbins, BINSIZE, 3,3,nsep).sum(axis=1)   # (nbins,3,3,nsep)
tot = binsum.sum(axis=0)
full = effmass_from_D(tot / N)
jk = np.zeros((nbins,3,3,len(tis)))
for b in range(nbins):
    jk[b] = effmass_from_D((tot - binsum[b]) / (N - BINSIZE))
mean = jk.mean(axis=0)
err = np.sqrt((nbins-1)/nbins * np.sum((jk-mean)**2, axis=0))

# write: t  then per (a,ilm) mean err
with open("shape_diag_effmass_Nf4g2L1_claude.dat","w") as f:
    hdr = "# t"
    for a in range(3):
        for ilm in range(3):
            hdr += f"  {SHAPES[a]}_{MCH[ilm]}_m {SHAPES[a]}_{MCH[ilm]}_e"
    f.write(hdr+"\n")
    for n,ti in enumerate(tis):
        row = f"{ti*AT:.4f}"
        for a in range(3):
            for ilm in range(3):
                row += f" {mean[a,ilm,n]:.6f} {err[a,ilm,n]:.6f}"
        f.write(row+"\n")
print("# wrote shape_diag_effmass_Nf4g2L1_claude.dat")
# quick spread check at small t
for n,ti in [(0,tis[0]),(2,tis[2])]:
    vals = mean[:,:,n].ravel()
    print(f"# t={ti*AT:.1f}: 9 (shape,m) values min={np.nanmin(vals):.4f} max={np.nanmax(vals):.4f} spread={np.nanmax(vals)-np.nanmin(vals):.4f}")
