#!/usr/bin/env python3
# eig_dov_spectrum_plot_claude.py
# _claude: scatter the dense complex spectrum of D_ov on the complex plane; overlay the GW circle |z-1|=1.
# Near-zero modes cluster near z=0 (left edge of the circle).  Usage: python3 THIS <eig_dov_spectrum_*.dat>

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "eig_dov_spectrum_L1_nt128_claude.dat"

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

n_near = int(np.sum(absz < 0.2))
imin = int(np.argmin(absz))

fig, ax = plt.subplots(figsize=(6.6, 6.4))

# GW circle |z-1| = 1
th = np.linspace(0.0, 2.0 * np.pi, 400)
ax.plot(1.0 + np.cos(th), np.sin(th), color="0.6", lw=1.0, ls="--", label=r"GW circle $|z-1|=1$")

ax.scatter(re, im, s=6, color="C0", alpha=0.5, label=r"eig$(D_{ov})$, $N=%d$" % len(re))
ax.scatter([re[imin]], [im[imin]], s=60, facecolors="none", edgecolors="C3", lw=1.8,
           label=r"smallest $|z|=%.4f$" % absz[imin])

ax.axvline(0.0, color="0.8", lw=0.8)
ax.axhline(0.0, color="0.8", lw=0.8)
ax.set_xlabel(r"$\mathrm{Re}\,z$")
ax.set_ylabel(r"$\mathrm{Im}\,z$")
ax.set_title(r"$D_{ov}$ complex spectrum (%s)" % path.replace("_claude.dat", ""))
ax.set_aspect("equal", "box")
ax.grid(True, alpha=0.25)
ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

fig.tight_layout()
out = "eig_dov_spectrum_claude.png"
fig.savefig(out, dpi=130)
print("wrote %s" % out)
print("N = %d  eigenvalues" % len(re))
print("smallest |z| = %.5f   (Re=%.5f, Im=%.5f)" % (absz[imin], re[imin], im[imin]))
print("largest  |z| = %.5f" % absz.max())
print("# modes with |z| < 0.2 : %d" % n_near)
