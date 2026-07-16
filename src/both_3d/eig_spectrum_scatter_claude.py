#!/usr/bin/env python3
# eig_spectrum_scatter_claude.py
# _claude: generic complex-spectrum scatter (Re vs Im) for any eig_*_spectrum_*.dat, with an optional
# vertical line marking the overlap wall Re = M (to see whether M sits in a clean gap between the physical
# branch near Re=0 and the Wilson doublers).  Usage:
#   python3 THIS <dat> <out.png> "<title>" [wall_M]

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else "eig_spectrum_scatter_claude.png"
title = sys.argv[3] if len(sys.argv) > 3 else path
wall = float(sys.argv[4]) if len(sys.argv) > 4 else None
zmax = float(sys.argv[5]) if len(sys.argv) > 5 else None   # drop |z|>zmax (spurious 1/theta Arnoldi outliers)

re = []
im = []
with open(path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        re.append(float(p[1]))
        im.append(float(p[2]))
re = np.array(re)
im = np.array(im)
absz = np.sqrt(re * re + im * im)
if zmax is not None:
    mask = absz < zmax
    re = re[mask]
    im = im[mask]
    absz = absz[mask]
    print("# kept %d of %d points with |z| < %.3g" % (len(re), len(mask), zmax))
imin = int(np.argmin(absz))

fig, ax = plt.subplots(figsize=(7.2, 6.2))
ax.scatter(re, im, s=8, color="C0", alpha=0.6, label=r"eig, $N=%d$" % len(re))
ax.scatter([re[imin]], [im[imin]], s=70, facecolors="none", edgecolors="C3", lw=1.8,
           label=r"smallest $|\lambda|=%.3f$" % absz[imin])
ax.axvline(0.0, color="0.8", lw=0.8)
ax.axhline(0.0, color="0.8", lw=0.8)
if wall is not None:
    ax.axvline(wall, color="C2", ls="--", lw=1.4, label=r"wall $\mathrm{Re}=M=%.2f$" % wall)

ax.set_xlabel(r"$\mathrm{Re}\,\lambda$")
ax.set_ylabel(r"$\mathrm{Im}\,\lambda$")
ax.set_title(title)
ax.set_xlim(-1.0, 4.0)          # FIXED range for the bare D_W low modes (origin, wall M=1, physical branch)
ax.set_ylim(-2.5, 2.5)
ax.set_aspect("equal", "box")
ax.grid(True, alpha=0.25)
ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote %s" % out)
print("N=%d  Re in [%.4f, %.4f]  smallest |lambda|=%.5f (Re=%.4f Im=%.4f)"
      % (len(re), re.min(), re.max(), absz[imin], re[imin], im[imin]))
