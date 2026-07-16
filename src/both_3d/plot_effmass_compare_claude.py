#!/usr/bin/env python3
# Effective-mass Delta_eff(t) curves for the shape-basis glueball GEVP (Nf2 g8), overlaying the
# operator bases 1/4/5/7 shapes so the plateau quality is directly visible. One figure per (channel,L).
#   F   (l=1)  : ground effmass, UNIT mean over the 3 m-blocks (degenerate triplet), jackknife error.
#   F^2 (0++)  : state 0 (heavier of nops2=2 per-m GEVP), jackknife error.
# Horizontal band = diagonal-chi2 const fit (F [0.2,0.8], F^2 [0.2,0.6]). Anchors: F->sqrt2, F^2->2sqrt2.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"
bases = ["tri", "4sh", "5sh", "7sh"]
# basis: distinct color AND marker
bstyle = {"tri": ("0.45", "x", "1 (triangle)"),
          "4sh": ("#d62728", "o", "4 (old prod)"),
          "5sh": ("#1f77b4", "s", "5 (+star)"),
          "7sh": ("#2ca02c", "^", "7 (+trio,five-six)")}
bdodge = {"tri": -0.03, "4sh": -0.01, "5sh": 0.01, "7sh": 0.03}

def load(fn):
    h = open(fn).readline().split()
    nb, nt, ns, at = int(h[1]), int(h[2]), int(h[3]), float(h[4])
    a = np.loadtxt(fn).reshape(nb, nt, ns)
    return a, at

def effmass_curve(jk, chan):
    # returns t (physical), mean effmass, jackknife error, per timeslice
    nb, nt, ns = jk.shape
    if chan == "F":
        # unit mean over the m-blocks per bin, then jackknife across bins
        c = jk.mean(axis=2)          # (nb, nt)
    else:
        c = jk[:, :, 0]              # state 0 (0++)
    mean = c.mean(0)
    err = np.sqrt((nb-1)/nb*np.sum((c-mean)**2, 0))
    return mean, err

def fit_band(jk, chan, tlo, thi, at):
    nb, nt, ns = jk.shape
    ti = [t for t in range(nt) if tlo-1e-9 <= (t+1)*at <= thi+1e-9]
    if chan == "F":
        # inv-var m-average matches the canonical fit
        Mbs = []
        sigs = []
        for s in range(ns):
            d = jk[:, ti, s]
            mn = d.mean(0)
            var = (nb-1)/nb*np.sum((d-mn)**2, 0)
            w = 1.0/var
            w /= w.sum()
            Mbs.append(d @ w)
            sigs.append(np.sqrt((nb-1)/nb*np.sum((Mbs[-1]-Mbs[-1].mean())**2)))
        Mbs = np.array(Mbs)
        u = np.array([1.0/s**2 for s in sigs])
        u /= u.sum()
        Mc = (u[:, None]*Mbs).sum(0)
        return Mc.mean(), np.sqrt((nb-1)/nb*np.sum((Mc-Mc.mean())**2))
    d = jk[:, ti, 0]
    mn = d.mean(0)
    var = (nb-1)/nb*np.sum((d-mn)**2, 0)
    w = 1.0/var
    w /= w.sum()
    Mb = d @ w
    return Mb.mean(), np.sqrt((nb-1)/nb*np.sum((Mb-Mb.mean())**2))

panels = [("F", 1, (0.2, 0.8), np.sqrt(2), r"free cont. $\sqrt{2}$", (1.0, 1.7), 8),
          ("F", 2, (0.2, 0.8), np.sqrt(2), r"free cont. $\sqrt{2}$", (0.8, 1.8), 8),
          ("F2", 1, (0.2, 0.6), 2*np.sqrt(2), r"two-photon $2\sqrt{2}$", (1.8, 3.6), 7),
          ("F2", 2, (0.2, 0.6), 2*np.sqrt(2), r"two-photon $2\sqrt{2}$", (0.0, 7.0), 7)]

for chan, L, rng, anchor, alab, ylim, ntshow in panels:
    plt.figure(figsize=(7.4, 5.0))
    ax = plt.gca()
    ax.axhline(anchor, ls="--", color="0.35", lw=1.1)
    ax.text(0.22, anchor, alab, ha="left", va="bottom", color="0.35", fontsize=9)
    ax.axvspan(rng[0], rng[1], color="0.85", alpha=0.5, zorder=0)
    for b in bases:
        if chan == "F2" and b == "tri":
            continue   # triangle cannot separate the vacuum constant for 0++
        fn = f"{SC}/cmp_{chan}_Nf2_g8.000000_L{L}_{b}.dat"
        jk, at = load(fn)
        mean, err = effmass_curve(jk, chan)
        c, mk, lab = bstyle[b]
        nt = min(ntshow, len(mean))
        tt = np.array([(i+1)*at for i in range(nt)]) + bdodge[b]
        ax.errorbar(tt, mean[:nt], yerr=err[:nt], color=c, marker=mk, ms=6, capsize=2,
                    lw=1.2, ls="-", label=lab, zorder=3)
        M, e = fit_band(jk, chan, rng[0], rng[1], at)
        ax.axhspan(M-e, M+e, xmin=(rng[0]-0.0)/(ntshow*at), xmax=(rng[1])/(ntshow*at),
                   color=c, alpha=0.12, zorder=1)
    ax.set_xlim(0.05, ntshow*at+0.05)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$t$  [$a_t$ units, $t = n\,a_t$]")
    ax.set_ylabel(r"$\Delta_\mathrm{eff}(t)$  [1/$a_t$]")
    ttl = "linear-$F$ $\\ell=1$" if chan == "F" else "$F^2$ 0++"
    ax.set_title(f"Effective mass, Nf2 $g^2$=8, L={L}: {ttl}  (shaded = fit range)")
    ax.legend(loc="upper right", frameon=True, fontsize=9)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    out = f"glue_{chan}_effmass_L{L}_claude.png"
    plt.savefig(out, dpi=140)
    plt.close()
    print("wrote", out)
