# plot_f2_spectra_L2_claude.py
# F^2 0++ GEVP effmass vs t at L2 for g^2 = 1,2,3 (Nf2/4/6), to see the large stat error at g^2=3.
# 0++ = GEVP state 0 (cols 3,4) of the F^2-v2 nops2=2 dump.  Fit window t=[0.2,0.4] shaded.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
GS = [1.0, 2.0, 3.0]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
TLO, THI = 0.2, 0.4

fig, axs = plt.subplots(1, len(GS), figsize=(6.5 * len(GS), 5), sharey=True)
for ax, g in zip(np.atleast_1d(axs).ravel(), GS):
    for nf in NFS:
        dat = "gevp_f2_Nf%d_g%.6f_L2_claude.dat" % (nf, g)
        try:
            d = np.loadtxt(dat)
        except OSError:
            continue
        t = d[:, 0]
        m = d[:, 3]        # 0++ = state 0 mean
        e = d[:, 4]
        sel = t <= 1.2
        ax.errorbar(t[sel], m[sel], e[sel], marker="o", ms=4, capsize=2, lw=0.8,
                    color=nfcol[nf], label="Nf%d" % nf)
    ax.axvspan(TLO, THI, color="gray", alpha=0.12)
    ax.set_title(r"$F^2$ $0^{++}$  L2 $g^2$=%.1f" % g)
    ax.set_xlabel(r"physical $t = a_t n_t$")
    ax.set_ylabel(r"$m_\mathrm{eff}$ (physical)")
    ax.set_xlim(0, 1.2)
    ax.set_ylim(0, 6)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
fig.suptitle(r"$F^2$ $0^{++}$ GEVP effmass at L2, $g^2$=1,2,3 -- fit window t[0.2,0.4] shaded (only 2 pts; g^2=3 noisy)")
fig.tight_layout()
fig.savefig("figs/f2_spectra_L2_claude.png", dpi=150)
print("# wrote f2_spectra_L2_claude.png")

# print the fit-window points
for g in GS:
    for nf in NFS:
        try:
            d = np.loadtxt("gevp_f2_Nf%d_g%.6f_L2_claude.dat" % (nf, g))
        except OSError:
            continue
        pts = [(d[k, 0], d[k, 3], d[k, 4]) for k in range(len(d)) if TLO - 1e-9 <= d[k, 0] <= THI + 1e-9]
        s = "  ".join("t%.1f=%.3f(%.3f)" % (t, m, e) for t, m, e in pts)
        print("  L2 g%.1f Nf%d : %s" % (g, nf, s))
