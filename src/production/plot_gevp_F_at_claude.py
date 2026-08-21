# plot_gevp_F_at_claude.py
# F (l=1) GEVP effective-mass spectra, a_t=0.2 vs a_t=0.1 OVERLAID, paired ensembles (L1 g1.0, L2 g2.0).
# Plotted curve = the m-AVERAGE (inverse-variance weighted per-m grounds), same as plot_glue_gevp.
# fit band M+/-sig from fit_perm variance-avg over t=[0.2,1.0].  Both a_t on physical t = a_t*n_t.
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
PAIRS = [(1, 1.000000), (2, 2.000000)]
TLO, THI = 0.2, 1.0


def mavg_curve(d):
    # inverse-variance weighted mean over the per-m ground cols (means 3,5,..; errs 4,6,..)
    nst = (d.shape[1] - 3) // 2
    ms = d[:, [3 + 2 * s for s in range(nst)]]
    es = d[:, [4 + 2 * s for s in range(nst)]]
    w = 1.0 / np.maximum(es, 1e-12) ** 2
    m = (ms * w).sum(1) / w.sum(1)
    e = np.sqrt(1.0 / w.sum(1))
    return m, e


def fit(jk):
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(TLO), str(THI), "0,1,2"],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "variance-avg" in line:
            seg = line.split("M =")[-1]
            m = seg.strip().split()[0]
            try:
                return float(m.split("(")[0]), float(m.split("(")[1].rstrip(")"))
            except (ValueError, IndexError):
                return None, None
    return None, None


fig, axs = plt.subplots(1, len(PAIRS), figsize=(7 * len(PAIRS), 5.5), sharey=True)
for ax, (L, g) in zip(np.atleast_1d(axs).ravel(), PAIRS):
    for nf in NFS:
        for at, suf, mk, fill in [(0.2, "", "o", True), (0.1, "_at0.1", "s", False)]:
            tag = "Nf%d_g%.6f_L%d%s" % (nf, g, L, suf)
            dat = "gevp_F_%s_claude.dat" % tag
            jk = "gevp_F_%s_jk_claude.dat" % tag
            try:
                d = np.loadtxt(dat)
            except OSError:
                continue
            t = d[:, 0]
            m, e = mavg_curve(d)
            sel = t <= 2.0
            ax.errorbar(t[sel], m[sel], e[sel], marker=mk, ms=4, lw=0.8, capsize=2,
                        color=nfcol[nf], mfc=nfcol[nf] if fill else "none",
                        ls="-" if fill else "--", label="Nf%d at%.1f" % (nf, at))
            M, sig = fit(jk)
            if M is not None:
                ax.hlines(M, TLO, THI, color=nfcol[nf], lw=1.4, alpha=0.6,
                          ls="-" if fill else "--")
    ax.axvspan(TLO, THI, color="gray", alpha=0.08)
    ax.set_title(r"F ($\ell{=}1$)  L%d $g^2$=%.1f" % (L, g))
    ax.set_xlabel(r"physical $t=a_t n_t$")
    ax.set_xlim(0, 2.0)
    ax.set_ylim(0, 4)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
np.atleast_1d(axs).ravel()[0].set_ylabel(r"$m_\mathrm{eff}$ (physical)")
fig.suptitle(r"F ($\ell{=}1$) GEVP effmass: $a_t{=}0.2$ (filled) vs $a_t{=}0.1$ (open) -- F is $a_t$-stable")
fig.tight_layout()
fig.savefig("figs/gevp_F_at_spectra_claude.png", dpi=150)
print("# wrote gevp_F_at_spectra_claude.png")
