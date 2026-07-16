#!/usr/bin/env python3
# Effective-mass Delta_eff(t) for the L1 F^2 0++ GEVP (Nf2 g8, 7-shape, state0), overlaying the four
# rebase locations t0 == trebase = 0,1,2,3. Shows how the rebase timeslice shapes the plateau.
# Horizontal bands = diagonal-chi2 const fit over [0.2,0.6]. Anchor: 2 sqrt2 (two-photon).
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"
# trebase: distinct color AND marker
trstyle = {0: ("#d62728", "o"), 1: ("#1f77b4", "s"), 2: ("#2ca02c", "^"), 3: ("#9467bd", "D")}
trdodge = {0: -0.03, 1: -0.01, 2: 0.01, 3: 0.03}
rng = (0.2, 0.6)
ntshow = 7

def load(fn):
    h = open(fn).readline().split()
    nb, nt, ns, at = int(h[1]), int(h[2]), int(h[3]), float(h[4])
    return np.loadtxt(fn).reshape(nb, nt, ns), at

def fit0(jk, tlo, thi, at):
    nb = jk.shape[0]
    ti = [t for t in range(jk.shape[1]) if tlo-1e-9 <= (t+1)*at <= thi+1e-9]
    d = jk[:, ti, 0]
    mean = d.mean(0)
    var = (nb-1)/nb*np.sum((d-mean)**2, 0)
    w = (1.0/var)
    w /= w.sum()
    Mb = d @ w
    return Mb.mean(), np.sqrt((nb-1)/nb*np.sum((Mb-Mb.mean())**2))

plt.figure(figsize=(7.4, 5.0))
ax = plt.gca()
ax.axhline(2*np.sqrt(2), ls="--", color="0.35", lw=1.1)
ax.text(0.22, 2*np.sqrt(2)+0.02, r"two-photon $2\sqrt{2}$", ha="left", va="bottom", color="0.35", fontsize=9)
ax.axvspan(rng[0], rng[1], color="0.85", alpha=0.5, zorder=0)
for TR in (0, 1, 2, 3):
    jk, at = load(f"{SC}/tr_f2_L1_tr{TR}.dat")
    nb, nt, ns = jk.shape
    c = jk[:, :, 0]
    mean = c.mean(0)
    err = np.sqrt((nb-1)/nb*np.sum((c-mean)**2, 0))
    col, mk = trstyle[TR]
    n = min(ntshow, len(mean))
    tt = np.array([(i+1)*at for i in range(n)]) + trdodge[TR]
    ax.errorbar(tt, mean[:n], yerr=err[:n], color=col, marker=mk, ms=6, capsize=2,
                lw=1.2, ls="-", label=f"$t_0$={TR} (rebase at $t={TR*at:.1f}$)", zorder=3)
    M, e = fit0(jk, rng[0], rng[1], at)
    ax.axhspan(M-e, M+e, color=col, alpha=0.10, zorder=1)
ax.set_xlim(0.05, ntshow*at+0.05)
ax.set_ylim(1.8, 3.6)
ax.set_xlabel(r"$t$  [$a_t$ units, $t = n\,a_t$]")
ax.set_ylabel(r"$\Delta_\mathrm{eff}(t)$  [1/$a_t$]")
ax.set_title(r"Effective mass, Nf2 $g^2$=8, L=1: $F^2$ 0++ (7-shape) -- rebase-$t_0$ overlay")
ax.legend(loc="upper right", frameon=True, fontsize=9)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("glue_F2_trebase_effmass_L1_claude.png", dpi=140)
plt.close()
print("wrote glue_F2_trebase_effmass_L1_claude.png")
