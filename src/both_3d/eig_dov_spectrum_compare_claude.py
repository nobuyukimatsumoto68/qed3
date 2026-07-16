#!/usr/bin/env python3
# eig_dov_spectrum_compare_claude.py
# _claude: compare the dense complex D_ov spectra across couplings (gsq = 2, 4, 8) at L1, one panel each, on
# the GW circle |z-1|=1.  Shows whether the empty near-zero (left) arc fills in as the coupling weakens.
# Reads eig_dov_spectrum_L1_gsq{2,4,8}_nt128_claude.dat.

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

gsqs = [2, 4, 8]
cols = ["C0", "C2", "C1"]


def load(path):
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
    return re, im


fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.2), sharex=True, sharey=True)
th = np.linspace(0.0, 2.0 * np.pi, 400)

for ax, g, col in zip(axes, gsqs, cols):
    path = "eig_dov_spectrum_L1_gsq%d_nt128_claude.dat" % g
    try:
        re, im = load(path)
    except FileNotFoundError:
        ax.set_title(r"gsq=%d  (missing %s)" % (g, path))
        continue
    absz = np.sqrt(re * re + im * im)
    imin = int(np.argmin(absz))
    n_near = int(np.sum(absz < 0.2))
    ax.plot(1.0 + np.cos(th), np.sin(th), color="0.6", lw=1.0, ls="--")
    ax.scatter(re, im, s=6, color=col, alpha=0.5)
    ax.scatter([re[imin]], [im[imin]], s=70, facecolors="none", edgecolors="C3", lw=1.8)
    ax.axvline(0.0, color="0.85", lw=0.8)
    ax.axhline(0.0, color="0.85", lw=0.8)
    ax.set_title(r"gsq=%d :  min$|z|$=%.3f,  #$|z|{<}0.2$=%d" % (g, absz[imin], n_near))
    ax.set_xlabel(r"$\mathrm{Re}\,z$")
    ax.set_aspect("equal", "box")
    ax.grid(True, alpha=0.25)
    print("gsq=%d : N=%d  min|z|=%.5f (Re=%.4f Im=%.4f)  max|z|=%.4f  #(|z|<0.2)=%d"
          % (g, len(re), absz[imin], re[imin], im[imin], absz.max(), n_near))

axes[0].set_ylabel(r"$\mathrm{Im}\,z$")
fig.suptitle(r"$D_{ov}$ complex spectrum vs coupling (L1, thermalized, frozen window)")
fig.tight_layout()
out = "eig_dov_spectrum_compare_claude.png"
fig.savefig(out, dpi=130)
print("wrote %s" % out)
