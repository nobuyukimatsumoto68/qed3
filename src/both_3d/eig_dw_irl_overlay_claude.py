#!/usr/bin/env python3
# eig_dw_irl_overlay_claude.py
# _claude: validate the IRL+Rayleigh-Ritz low D_W spectrum against the DENSE reference (free L1).  Overlays
# the dense bare D_W eigenvalues (gray, = M_DW + M) and the IRL Ritz values (red rings, already bare) on the
# complex plane with the domain wall M.  If the red rings sit on the gray low-mode points, the method works.
# Usage: python3 THIS <dense_MDW.dat> <irl_bareDW.dat> <out.png> [M]

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dense_path = sys.argv[1]
irl_path = sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else "eig_dw_irl_overlay_claude.png"
M = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0


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
    return np.array(re), np.array(im)


# dense stores M_DW = D_W - M -> shift +M to the bare frame; IRL already stores bare D_W.
d_re, d_im = load(dense_path)
d_re = d_re + M
i_re, i_im = load(irl_path)
i_abs = np.sqrt(i_re * i_re + i_im * i_im)
iw = int(np.argmin(i_abs))

fig, ax = plt.subplots(figsize=(7.6, 6.6))
ax.scatter(d_re, d_im, s=10, color="0.6", alpha=0.6, label=r"dense $D_W$ ($N=%d$)" % len(d_re))
ax.scatter(i_re, i_im, s=55, facecolors="none", edgecolors="C3", lw=1.5,
           label=r"IRL Ritz ($|D_W|{<}2$, %d modes)" % len(i_re))
ax.axvline(M, color="C2", ls="--", lw=1.5, label=r"wall $M=%.1f$" % M)
ax.axvline(0.0, color="0.85", lw=0.8)
ax.axhline(0.0, color="0.85", lw=0.8)

ax.set_xlim(-0.5, 4.0)          # zoom on the low region (physical branch + first doublers)
ax.set_ylim(-4.0, 4.0)
ax.set_xlabel(r"$\mathrm{Re}\,D_W$")
ax.set_ylabel(r"$\mathrm{Im}\,D_W$")
ax.set_title(r"IRL low $D_W$ modes vs DENSE (free L1);  smallest $|D_W|_{IRL}=%.3f$" % i_abs[iw])
ax.grid(True, alpha=0.25)
ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote %s" % out)
print("dense: N=%d ; IRL: %d modes, smallest |D_W|=%.5f (Re=%.4f Im=%.4f)"
      % (len(d_re), len(i_re), i_abs[iw], i_re[iw], i_im[iw]))
