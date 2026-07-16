#!/usr/bin/env python3
# Side-by-side basis comparison of the shape-basis glueball masses (Nf2 g8), same configs, varying
# only the GEVP operator basis: 1 (triangle), 4 (old prod {tri,rect,fig8,three-tri}), 5 (+star),
# 7 (+trio,+five-six). Answers "did the new shapes improve the signal?". Series = L1 (filled) / L2 (open).
#   F   (l=1)  : inv-var m-average of per-m GEVP ground, diagonal-chi2 const fit [0.2,0.8]
#   F^2 (0++)  : state 0 (heavier of nops2=2), diagonal-chi2 const fit [0.2,0.6]
# Anchors: F l=1 -> sqrt2 ; F^2 0++ -> 2 sqrt2. Triangle-only F^2 CANNOT separate the vacuum
# constant from the 0++ (returns ~0), so it is shown as an annotated "fails" marker off-scale-clamped.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"
bases = ["tri", "4sh", "5sh", "7sh"]
blab = ["1\n(tri)", "4\n(old prod)", "5\n(+star)", "7\n(+trio,\nfive-six)"]

def load(fn):
    h = open(fn).readline().split()
    nb, nt, ns, at = int(h[1]), int(h[2]), int(h[3]), float(h[4])
    a = np.loadtxt(fn).reshape(nb, nt, ns)
    return a, at

def fit_state(jk, tlo, thi, at, s):
    nb = jk.shape[0]
    ti = [t for t in range(jk.shape[1]) if tlo-1e-9 <= (t+1)*at <= thi+1e-9]
    d = jk[:, ti, s]
    mean = d.mean(0)
    var = (nb-1)/nb*np.sum((d-mean)**2, 0)
    w = (1.0/var)
    w /= w.sum()
    Mb = d @ w
    M = Mb.mean()
    sig = np.sqrt((nb-1)/nb*np.sum((Mb-M)**2))
    return M, sig

def favg(jk, tlo, thi, at):
    nst = jk.shape[2]
    Mbs = []
    sigs = []
    nb = jk.shape[0]
    ti = [t for t in range(jk.shape[1]) if tlo-1e-9 <= (t+1)*at <= thi+1e-9]
    for s in range(nst):
        d = jk[:, ti, s]
        mean = d.mean(0)
        var = (nb-1)/nb*np.sum((d-mean)**2, 0)
        w = (1.0/var)
        w /= w.sum()
        Mbs.append(d @ w)
        sigs.append(np.sqrt((nb-1)/nb*np.sum((Mbs[-1]-Mbs[-1].mean())**2)))
    Mbs = np.array(Mbs)
    u = np.array([1.0/s**2 for s in sigs])
    u /= u.sum()
    Mc = (u[:, None]*Mbs).sum(0)
    return Mc.mean(), np.sqrt((nb-1)/nb*np.sum((Mc-Mc.mean())**2))

def get(chan, L, b, rng):
    fn = f"{SC}/cmp_{chan}_Nf2_g8.000000_L{L}_{b}.dat"
    jk, at = load(fn)
    if chan == "F":
        return favg(jk, rng[0], rng[1], at)
    return fit_state(jk, rng[0], rng[1], at, 0)

# L1 filled red circle, L2 open blue square (distinct color AND marker)
Lstyle = {1: ("#d62728", "o", True), 2: ("#1f77b4", "s", False)}
dodge = {1: -0.09, 2: 0.09}
x = np.arange(len(bases))

# ============================ Figure 1: F (l=1) ============================
plt.figure(figsize=(7.2, 5.0))
ax = plt.gca()
ax.axhline(np.sqrt(2), ls="--", color="0.35", lw=1.2)
ax.text(3.4, np.sqrt(2)+0.004, r"free cont. $\sqrt{2}$", ha="right", va="bottom", color="0.35", fontsize=9)
ax.axhline(1.33242, ls=":", color="0.55", lw=1.0)
ax.text(3.4, 1.33242-0.006, r"free lattice (L=1) det", ha="right", va="top", color="0.55", fontsize=9)
for L in (1, 2):
    c, mk, fill = Lstyle[L]
    for i, b in enumerate(bases):
        M, e = get("F", L, b, (0.2, 0.8))
        ax.errorbar(i+dodge[L], M, yerr=e, color=c, marker=mk, ms=8, capsize=3,
                    mfc=(c if fill else "white"), mec=c, lw=1.4, zorder=3)
        ax.annotate(f"{M:.3f}", (i+dodge[L], M+e), textcoords="offset points",
                    xytext=(0, 5), ha="center", va="bottom", fontsize=7.5, color=c)
    ax.plot([], [], color=c, marker=mk, ms=8, lw=1.4,
            mfc=(c if fill else "white"), mec=c, label=f"L={L}")
ax.set_xticks(x)
ax.set_xticklabels(blab)
ax.set_xlabel("GEVP basis (number of shapes)")
ax.set_ylabel(r"$M_F$  (linear $F$, $\ell=1$)  [1/$a_t$]")
ax.set_title(r"Basis comparison, Nf2 $g^2$=8: linear-$F$ $\ell=1$")
ax.set_xlim(-0.5, 3.5)
ax.legend(loc="lower right", frameon=True)
ax.grid(alpha=0.25, axis="y")
plt.tight_layout()
plt.savefig("glue_F_basiscmp_claude.png", dpi=140)
plt.close()

# ============================ Figure 2: F^2 (0++) ============================
plt.figure(figsize=(7.2, 5.0))
ax = plt.gca()
ax.axhline(2*np.sqrt(2), ls="--", color="0.35", lw=1.2)
ax.text(3.4, 2*np.sqrt(2)+0.02, r"free cont. two-photon $2\sqrt{2}$", ha="right", va="bottom", color="0.35", fontsize=9)
for L in (1, 2):
    c, mk, fill = Lstyle[L]
    for i, b in enumerate(bases):
        M, e = get("F2", L, b, (0.2, 0.6))
        if b == "tri":
            # triangle alone cannot separate the vacuum constant -> M ~ 0; mark as failure
            if L == 1:
                ax.annotate("1-shape: fails\n(vacuum const.\nnot separated)", (0, 3.15),
                            ha="center", va="center", fontsize=8, color="0.25",
                            bbox=dict(boxstyle="round", fc="#fff3f3", ec="0.5", alpha=0.95))
            continue
        ax.errorbar(i+dodge[L], M, yerr=e, color=c, marker=mk, ms=8, capsize=3,
                    mfc=(c if fill else "white"), mec=c, lw=1.4, zorder=3)
        ax.annotate(f"{M:.3f}", (i+dodge[L], M+e), textcoords="offset points",
                    xytext=(0, 5), ha="center", va="bottom", fontsize=7.5, color=c)
    ax.plot([], [], color=c, marker=mk, ms=8, lw=1.4,
            mfc=(c if fill else "white"), mec=c, label=f"L={L}")
ax.set_xticks(x)
ax.set_xticklabels(blab)
ax.set_xlabel("GEVP basis (number of shapes)")
ax.set_ylabel(r"$M_{F^2}$  (0++, $\ell=0$)  [1/$a_t$]")
ax.set_title(r"Basis comparison, Nf2 $g^2$=8: $F^2$ 0++")
ax.set_xlim(-0.5, 3.5)
ax.legend(loc="lower right", frameon=True)
ax.grid(alpha=0.25, axis="y")
plt.tight_layout()
plt.savefig("glue_F2_basiscmp_claude.png", dpi=140)
plt.close()
print("wrote glue_F_basiscmp_claude.png , glue_F2_basiscmp_claude.png")
