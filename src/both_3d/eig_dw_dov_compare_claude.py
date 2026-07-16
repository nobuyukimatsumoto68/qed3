#!/usr/bin/env python3
# eig_dw_dov_compare_claude.py
# _claude: side-by-side dense complex spectra of the bare Wilson D_W (left) and the overlap D_ov (right) for
# the SAME config.  D_ov is a UNITARY PROJECTOR of D_W: it keeps only each D_W eigenvalue's phase and maps it
# onto the GW circle |z-1|=1 -- so a small |D_W| eigenvalue (e.g. the 0.2 mode) does NOT give small |z|; only
# its phase matters.  Usage: python3 THIS <dw.dat> <dov.dat>

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dw_path = sys.argv[1] if len(sys.argv) > 1 else "eig_dw_spectrum_L1_nt128_claude.dat"
dov_path = sys.argv[2] if len(sys.argv) > 2 else "eig_dov_spectrum_L1_nt128_claude.dat"
out = sys.argv[3] if len(sys.argv) > 3 else "eig_dw_dov_compare_claude.png"
tag = sys.argv[4] if len(sys.argv) > 4 else ""


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


M_WALL = 1.0     # domain-wall height M.  Plotted data is M_DW = D_W - M; bare D_W = M_DW + M.

dw_re, dw_im = load(dw_path)
dov_re, dov_im = load(dov_path)
# recover the BARE Wilson D_W = M_DW + M (undo the -M baked into the stored M_DW).
dwb_re = dw_re + M_WALL
dw_abs = np.sqrt(dw_re * dw_re + dw_im * dw_im)          # |M_DW| = |A| = distance to the wall
dov_abs = np.sqrt(dov_re * dov_re + dov_im * dov_im)
iw = int(np.argmin(dw_abs))
iv = int(np.argmin(dov_abs))
n_light = int(np.sum(dwb_re < M_WALL))                    # modes on the light side Re(D_W) < M

fig, (axl, axr) = plt.subplots(1, 2, figsize=(12.0, 5.6))

# left: bare D_W eigenvalue cloud, with the domain wall M
axl.scatter(dwb_re, dw_im, s=6, color="C0", alpha=0.5)
axl.axvline(M_WALL, color="C2", ls="--", lw=1.6, label=r"domain wall $M=%.1f$" % M_WALL)
axl.scatter([dwb_re[iw]], [dw_im[iw]], s=70, facecolors="none", edgecolors="C3", lw=1.8,
            label=r"closest to wall $|A|=%.3f$ (#light=%d)" % (dw_abs[iw], n_light))
axl.axvline(0.0, color="0.8", lw=0.8)
axl.axhline(0.0, color="0.8", lw=0.8)
axl.set_title(r"bare $D_W$ eigenvalues;  Re$<M$ = light ($z\to0$),  Re$>M$ = doubler ($z\to2$)")
axl.set_xlabel(r"$\mathrm{Re}\,D_W$")
axl.set_ylabel(r"$\mathrm{Im}\,D_W$")
axl.set_aspect("equal", "box")
axl.grid(True, alpha=0.25)
axl.legend(loc="upper right", fontsize=8)

# right: D_ov on the GW circle
th = np.linspace(0.0, 2.0 * np.pi, 400)
axr.plot(1.0 + np.cos(th), np.sin(th), color="0.6", lw=1.0, ls="--", label=r"GW circle $|z-1|=1$")
axr.scatter(dov_re, dov_im, s=6, color="C1", alpha=0.5)
axr.scatter([dov_re[iv]], [dov_im[iv]], s=70, facecolors="none", edgecolors="C3", lw=1.8,
            label=r"smallest $|z|=%.3f$" % dov_abs[iv])
axr.axvline(0.0, color="0.8", lw=0.8)
axr.axhline(0.0, color="0.8", lw=0.8)
axr.set_title(r"$D_{ov}$ eigenvalues (unitary projection of $D_W$)")
axr.set_xlabel(r"$\mathrm{Re}\,z$")
axr.set_ylabel(r"$\mathrm{Im}\,z$")
axr.set_aspect("equal", "box")
axr.grid(True, alpha=0.25)
axr.legend(loc="upper right", fontsize=8)

fig.suptitle(r"$D_W$ vs $D_{ov}$ %s: overlap keeps the PHASE, discards $|\lambda|$" % tag)
fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote %s" % out)
print("D_W  : N=%d  smallest |lambda|=%.5f  largest |lambda|=%.5f" % (len(dw_re), dw_abs[iw], dw_abs.max()))
print("D_ov : N=%d  smallest |z|=%.5f       largest |z|=%.5f" % (len(dov_re), dov_abs[iv], dov_abs.max()))
