#!/usr/bin/env python3
# Trend plots of the L=1 shape-basis glueball masses vs gsq, one series per Nf.
#   F   (l=1)  : inv-var m-average of the per-m GEVP ground, diagonal-chi2 const fit over [0.2,0.8]
#   F^2 (0++)  : state 0 (heavier of nops2=2), diagonal-chi2 const fit over [0.2,0.6]
# Free-continuum anchors: F l=1 -> sqrt2 ; F^2 0++ -> 2 sqrt2 (two-photon). L=1 lattice free det -> 1.33242.
# Reliable = gsq2, gsq8 (3k-6k cfg, filled markers); thin = gsq1,4,12 (319 cfg, open markers).
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"

def load(fn):
    h = open(fn).readline().split()
    nbins, ntpts, nst, at = int(h[1]), int(h[2]), int(h[3]), float(h[4])
    a = np.loadtxt(fn).reshape(nbins, ntpts, nst)
    return a, at

def fit_state(jk, tlo, thi, at, s):
    nb = jk.shape[0]
    ti = [t for t in range(jk.shape[1]) if tlo-1e-9 <= (t+1)*at <= thi+1e-9]
    d = jk[:, ti, s]
    mean = d.mean(0)
    dev = d - mean
    var = (nb-1)/nb * np.sum(dev*dev, 0)
    w = (1.0/var)
    w /= w.sum()
    Mb = d @ w
    M = Mb.mean()
    sig = np.sqrt((nb-1)/nb * np.sum((Mb-M)**2))
    return Mb, M, sig

def favg(jk, tlo, thi, at, states):
    Mbs = []
    sigs = []
    for s in states:
        Mb, M, sig = fit_state(jk, tlo, thi, at, s)
        Mbs.append(Mb)
        sigs.append(sig)
    Mbs = np.array(Mbs)
    u = np.array([1.0/s**2 for s in sigs])
    u /= u.sum()
    Mc = (u[:, None]*Mbs).sum(0)
    nb = jk.shape[0]
    M = Mc.mean()
    err = np.sqrt((nb-1)/nb * np.sum((Mc-M)**2))
    return M, err

gsqs = [("1.000000", 1), ("2.000000", 2), ("4.000000", 4), ("8.000000", 8), ("12.000000", 12)]
reliable = {2, 8}
# Nf series: distinct color AND marker (color-blind safe)
style = {2: ("#d62728", "o"), 4: ("#1f77b4", "s"), 6: ("#2ca02c", "^")}
# horizontal dodge per Nf (multiplicative on the log2 axis) so overlapping points/labels separate
dodge = {2: 1.0/1.09, 4: 1.0, 6: 1.09}

# ---- gather ----
res = {}   # res[(chan, NF)] = list of (gsq, M, err, reliable)
for chan in ("F", "F2"):
    for NF in (2, 4, 6):
        rows = []
        for g, gi in gsqs:
            tag = f"Nf{NF}_g{g}_L1"
            if chan == "F":
                jk, at = load(f"{SC}/jk_msm_{tag}.dat")
                M, e = favg(jk, 0.2, 0.8, at, [0, 1, 2])
            else:
                jk, at = load(f"{SC}/jk_f2b_{tag}.dat")
                _, M, e = fit_state(jk, 0.2, 0.6, at, 0)
            rows.append((gi, M, e, gi in reliable))
        res[(chan, NF)] = rows

# ============================ Figure 1: F (l=1) ============================
plt.figure(figsize=(7.2, 5.0))
ax = plt.gca()
ax.axhline(np.sqrt(2), ls="--", color="0.35", lw=1.2)
ax.text(0.95, np.sqrt(2)+0.004, r"free cont. $\sqrt{2}=1.4142$", ha="left", va="bottom", color="0.35", fontsize=9)
ax.axhline(1.33242, ls=":", color="0.55", lw=1.0)
ax.text(0.95, 1.33242-0.010, r"free lattice (L=1) det $=1.3324$", ha="left", va="top", color="0.55", fontsize=9)
for NF in (2, 4, 6):
    c, mk = style[NF]
    rows = res[("F", NF)]
    for x, y, e, f in rows:
        xp = x*dodge[NF]
        ax.errorbar(xp, y, yerr=e, color=c, marker=mk, ms=8, capsize=3,
                    mfc=(c if f else "white"), mec=c, lw=1.4, zorder=3)
        ax.annotate(f"{y:.3f}", (xp, y+e), textcoords="offset points", xytext=(0, 5),
                    ha="center", va="bottom", fontsize=7.5, color=c)
    ax.plot([], [], color=c, marker=mk, ms=8, lw=1.4, label=f"Nf={NF}")
ax.set_xscale("log", base=2)
ax.set_xticks([1, 2, 4, 8, 12])
ax.set_xticklabels(["1", "2", "4", "8", "12"])
ax.set_xlabel(r"$g^2$")
ax.set_ylabel(r"$M_F$  (linear $F$, $\ell=1$)  [1/$a_t$]")
ax.set_title(r"$L=1$ glueball: linear-$F$ $\ell=1$ mass vs $g^2$ (filled = reliable gsq2,8)")
ax.set_xlim(0.85, 14)
ax.legend(loc="lower left", frameon=True)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("glue_F_trend_L1_claude.png", dpi=140)
plt.close()

# ============================ Figure 2: F^2 (0++) ============================
plt.figure(figsize=(7.2, 5.0))
ax = plt.gca()
ax.axhline(2*np.sqrt(2), ls="--", color="0.35", lw=1.2)
ax.text(0.95, 2*np.sqrt(2)+0.015, r"free cont. two-photon $2\sqrt{2}=2.8284$", ha="left", va="bottom", color="0.35", fontsize=9)
for NF in (2, 4, 6):
    c, mk = style[NF]
    rows = res[("F2", NF)]
    for x, y, e, f in rows:
        xp = x*dodge[NF]
        ax.errorbar(xp, y, yerr=e, color=c, marker=mk, ms=8, capsize=3,
                    mfc=(c if f else "white"), mec=c, lw=1.4, zorder=3)
        ax.annotate(f"{y:.3f}", (xp, y+e), textcoords="offset points", xytext=(0, 5),
                    ha="center", va="bottom", fontsize=7.5, color=c)
    ax.plot([], [], color=c, marker=mk, ms=8, lw=1.4, label=f"Nf={NF}")
ax.set_xscale("log", base=2)
ax.set_xticks([1, 2, 4, 8, 12])
ax.set_xticklabels(["1", "2", "4", "8", "12"])
ax.set_xlabel(r"$g^2$")
ax.set_ylabel(r"$M_{F^2}$  (0++, $\ell=0$)  [1/$a_t$]")
ax.set_title(r"$L=1$ glueball: $F^2$ 0++ mass vs $g^2$ (filled = reliable gsq2,8)")
ax.set_xlim(0.85, 14)
ax.legend(loc="lower left", frameon=True)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("glue_F2_trend_L1_claude.png", dpi=140)
plt.close()
print("wrote glue_F_trend_L1_claude.png , glue_F2_trend_L1_claude.png")
