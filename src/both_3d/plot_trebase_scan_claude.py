#!/usr/bin/env python3
# Rebase-location (t0 == trebase) scan for the L1 F^2 0++ GEVP (Nf2 g8, 7-shape, state0).
# The fixed-t0 rebase uses (matCA[trebase], matCA[trebase+dt]); scanning trebase=0,1,2,3 tests the
# sensitivity of the extracted 0++ mass to the rebase timeslice. Fit = diagonal-chi2 const [0.2,0.6].
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"

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
    return M, sig, chi2, len(ti)-1

trs = [0, 1, 2, 3]
Ms = []
es = []
cs = []
for TR in trs:
    jk, at = load(f"{SC}/tr_f2_L1_tr{TR}.dat")
    M, e, c, d = fit0(jk, 0.2, 0.6, at)
    Ms.append(M)
    es.append(e)
    cs.append(c/max(d, 1))

plt.figure(figsize=(7.0, 5.0))
ax = plt.gca()
ax.axhline(2*np.sqrt(2), ls="--", color="0.35", lw=1.2)
ax.text(3.35, 2*np.sqrt(2)+0.006, r"free cont. two-photon $2\sqrt{2}=2.8284$", ha="right", va="bottom", color="0.35", fontsize=9)
c = "#9467bd"
ax.errorbar(trs, Ms, yerr=es, color=c, marker="D", ms=9, capsize=4, lw=1.5, ls="-", zorder=3)
for TR, M, e, ch in zip(trs, Ms, es, cs):
    ax.annotate(f"{M:.3f}\n$\\chi^2$/dof={ch:.2f}", (TR, M+e), textcoords="offset points",
                xytext=(0, 6), ha="center", va="bottom", fontsize=8, color=c)
ax.set_xticks(trs)
ax.set_xlabel(r"rebase location $t_0$ (= trebase, timeslice index; $t_0 a_t = 0,\,0.2,\,0.4,\,0.6$)")
ax.set_ylabel(r"$M_{F^2}$  (0++, $\ell=0$)  [1/$a_t$]")
ax.set_title(r"Rebase-$t_0$ scan, Nf2 $g^2$=8, L=1: $F^2$ 0++ (7-shape, fit [0.2,0.6])")
ax.set_xlim(-0.4, 3.4)
ax.set_ylim(2.55, 3.15)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("glue_F2_trebase_L1_claude.png", dpi=140)
plt.close()
print("wrote glue_F2_trebase_L1_claude.png")
