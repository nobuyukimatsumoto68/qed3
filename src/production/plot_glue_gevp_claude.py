# plot_glue_gevp_claude.py
# GEVP effective-mass spectra for the massless redo ensembles (per-m strategy, kmin=20).
# Two figures (F l=1 ground, F^2 0++), each panelled by L (1,2), curves = ensembles (Nf x gsq).
# Effmass points from gevp_<chan>_<tag>_claude.dat (col1=t, col2=ground mean, col3=err);
# fit band (M +/- sig) from fit_perm on the jk dump over the channel's t-window.
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NFS = [2, 4, 6]
GS = {1: [0.5, 1.0, 1.5], 2: [1.0, 2.0, 3.0], 3: [1.5, 3.0, 4.5], 4: [2.0, 4.0, 6.0]}
nfcol = {2: "tab:red", 4: "tab:blue", 6: "tab:green"}
gls = {0.5: "-", 1.0: "--", 1.5: "-", 2.0: "--", 3.0: ":", 4.0: "--", 4.5: ":", 6.0: ":"}
# channel: (dat prefix, fit t-window, fit states, plotted-state col pair, title)
# .dat cols: 0=t, 1/2=ground(=last state=lightest), 3/4=s0, 5/6=s1, 7/8=s2, 9/10=s3, ...
# F l=1 ground = the m-averaged plateau, use the ground col (1,2).
# F^2 0++ INCLUDES F^4 (NM 2026-08-18): prefix glue_f2_v2_shapes, nops2=4 -> 4 states.  With vacsub=0
# the near-zero VACUUM mode is the lightest (ground col / state 3), so the physical 0++ is STATE 2 =
# cols (7,8).  (Old 7-shape F^2 was nops2=2 -> 0++ = s0 = cols (3,4).)
CH = {
    "F":   ("gevp_F",   (0.2, 1.0), "0,1,2",     (1, 2), r"$F$ ($\ell=1$) ground"),
    "Fl2": ("gevp_Fl2", (0.2, 0.6), "0,1,2,3,4", (1, 2), r"$F$ ($\ell=2$, H) ground"),
    "f2":  ("gevp_f2",  (0.2, 0.4), "0",         (3, 4), r"$F^2$/$F^4$ ($0^{++}$)"),
}


def fit(jk, tlo, thi, states):
    out = subprocess.run(["python3", "fit_perm_claude.py", jk, str(tlo), str(thi), states],
                         capture_output=True, text=True).stdout
    single = ("," not in states)   # single-state fit (e.g. F^2 0++ = "2"): match "state <s>:"
    for line in out.splitlines():
        if "variance-avg" in line or (single and ("state %s:" % states) in line):
            # parse M(err)
            seg = line.split("M =")[-1] if "==>" in line else line.split("M=")[-1]
            m = seg.strip().split()[0]
            val = float(m.split("(")[0])
            err = float(m.split("(")[1].rstrip(")"))
            return val, err
    return None, None


def plotted_curve(d, chan, cm, ce):
    # per-m channels (F, Fl2): the plotted curve must MATCH what is fitted -- the m-AVERAGE over the
    # per-m ground states, NOT the single-state "ground" column (col 1/2 = lightest m-block only).
    # inverse-variance weighted mean per t over the per-m state cols (means=3,5,..; errs=4,6,..).
    if chan in ("F", "Fl2"):
        nst = (d.shape[1] - 3) // 2
        ms = d[:, [3 + 2 * s for s in range(nst)]]
        es = d[:, [4 + 2 * s for s in range(nst)]]
        w = 1.0 / np.maximum(es, 1e-12) ** 2
        m = (ms * w).sum(1) / w.sum(1)
        e = np.sqrt(1.0 / w.sum(1))
        return m, e
    return d[:, cm], d[:, ce]


for chan, (pref, (tlo, thi), states, (cm, ce), title) in CH.items():
    # F l=1 and F l=2: all 4 L.  F^2/F^4 0++: L1,L2 ONLY (NM 2026-08-18 -- L3/L4 0++ unresolved at
    # these stats; 10/36 returned nan, all strong-coupling L3/L4).
    Ls = [1, 2] if chan == "f2" else [1, 2, 3]   # L=4 dropped for F (NM 2026-08-18)
    fig, axs = plt.subplots(1, len(Ls), figsize=(6 * len(Ls), 5), sharey=True)
    for ax, L in zip(np.atleast_1d(axs), Ls):
        for nf in NFS:
            for g in GS[L]:
                tag = "Nf%d_g%.6f_L%d" % (nf, g, L)
                dat = "%s_%s_claude.dat" % (pref, tag)
                jk = "%s_%s_jk_claude.dat" % (pref, tag)
                try:
                    d = np.loadtxt(dat)
                except OSError:
                    continue
                t = d[:, 0]
                m, e = plotted_curve(d, chan, cm, ce)
                sel = t <= 1.6
                lab = "Nf%d g%.1f" % (nf, g)
                ax.errorbar(t[sel], m[sel], e[sel], marker="o", ms=3, lw=0.8,
                            color=nfcol[nf], ls=gls[g], capsize=2, label=lab)
                M, sig = fit(jk, tlo, thi, states)
                if M is not None:
                    ax.hlines(M, tlo, thi, color=nfcol[nf], lw=1.5, alpha=0.7)
        ax.axvspan(tlo, thi, color="gray", alpha=0.08)
        ax.set_title("L%d  %s" % (L, title))
        ax.set_xlabel("t")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6, ncol=2)
        ax.set_ylim(0, 5)
    np.atleast_1d(axs)[0].set_ylabel(r"$a_t\, m_{\rm eff}$")
    fig.tight_layout()
    fig.savefig("figs/glue_gevp_%s_spectra_claude.png" % chan, dpi=150)
    print("# wrote glue_gevp_%s_spectra_claude.png" % chan)
