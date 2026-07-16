#!/usr/bin/env python3
# L2 F^2 0++ (Nf2 g8, t0=1, nops2=2, state0): effect of dropping fig8 and trio together, vs the full
# 7-op basis and the single drops. Dropping fig8(2)+trio(5) -> 5-op basis {tri,rect,three-tri,star,
# five-six}. Bands = diagonal-chi2 const fit [0.2,0.6]. Anchor 2 sqrt2.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"
variants = ["full", "drop_fig8_trio", "drop_fig8_trio_star"]
vlabel = {"full": "full (7 ops)", "drop_fig8_trio": "drop fig8+trio (5 ops)",
          "drop_fig8_trio_star": "drop fig8+trio+star (4 ops)"}
vstyle = {"full": ("k", "o"), "drop_fig8_trio": ("#2ca02c", "D"),
          "drop_fig8_trio_star": ("#d62728", "s")}
vdodge = {"full": -0.015, "drop_fig8_trio": 0.0, "drop_fig8_trio_star": 0.015}
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
    M = Mb.mean()
    sig = np.sqrt((nb-1)/nb*np.sum((Mb-M)**2))
    chi2 = float(np.sum((mean-M)**2/var))
    return M, sig, chi2/max(len(ti)-1, 1)

plt.figure(figsize=(7.6, 5.2))
ax = plt.gca()
ax.axhline(2*np.sqrt(2), ls="--", color="0.35", lw=1.1)
ax.text(0.22, 2*np.sqrt(2)+0.04, r"two-photon $2\sqrt{2}$", ha="left", va="bottom", color="0.35", fontsize=9)
ax.axvspan(rng[0], rng[1], color="0.85", alpha=0.5, zorder=0)
for v in variants:
    jk, at = load(f"{SC}/dg_f2_L2_{v}.dat")
    nb = jk.shape[0]
    c = jk[:, :, 0]
    mean = c.mean(0)
    err = np.sqrt((nb-1)/nb*np.sum((c-mean)**2, 0))
    col, mk = vstyle[v]
    n = min(ntshow, jk.shape[1])
    tt = np.array([(i+1)*at for i in range(n)]) + vdodge[v]
    M, e, ch = fit0(jk, rng[0], rng[1], at)
    ax.errorbar(tt, mean[:n], yerr=err[:n], color=col, marker=mk, ms=6, capsize=2, lw=1.8,
                ls="-", label=f"{vlabel[v]}:  M={M:.3f}({e:.3f}), $\\chi^2$/dof={ch:.2f}",
                zorder=3)
    ax.axhspan(M-e, M+e, color=col, alpha=0.10, zorder=1)
ax.set_xlim(0.05, ntshow*at+0.05)
ax.set_ylim(0.0, 7.0)
ax.set_xlabel(r"$t$  [$a_t$ units, $t = n\,a_t$]")
ax.set_ylabel(r"$\Delta_\mathrm{eff}(t)$  [1/$a_t$]")
ax.set_title(r"Effective mass, Nf2 $g^2$=8, L=2: $F^2$ 0++ -- progressive drop (variational bias up)")
ax.legend(loc="upper right", frameon=True, fontsize=8.5)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("glue_F2_dropfig8trio_effmass_L2_claude.png", dpi=140)
plt.close()
print("wrote glue_F2_dropfig8trio_effmass_L2_claude.png")
