#!/usr/bin/env python3
# nops2 (rebase-subspace size) scan for L1 F^2 0++ (Nf2 g8, 7-shape, trebase t0=1). Enlarging the
# kept rebase subspace nops2=2,3,4 is a variational check on the lightest 0++. In gevp_group the
# states come out heaviest->lightest, so the vacuum constant is the LAST state and the lightest 0++
# is the second-to-last (index nops2-2). Overlay its Delta_eff(t); nops2=3,4 also resolve a heavier
# excited 0++ (~5), plotted thin/dashed. Bands = diagonal-chi2 const fit [0.2,0.6]. Anchor 2 sqrt2.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"
nstyle = {2: ("#d62728", "o"), 3: ("#1f77b4", "s"), 4: ("#2ca02c", "^")}
ndodge = {2: -0.025, 3: 0.0, 4: 0.025}
rng = (0.2, 0.6)
ntshow = 7

def load(fn):
    h = open(fn).readline().split()
    nb, nt, ns, at = int(h[1]), int(h[2]), int(h[3]), float(h[4])
    return np.loadtxt(fn).reshape(nb, nt, ns), at

def curve(jk, s):
    nb = jk.shape[0]
    c = jk[:, :, s]
    mean = c.mean(0)
    err = np.sqrt((nb-1)/nb*np.sum((c-mean)**2, 0))
    return mean, err

def fit(jk, s, tlo, thi, at):
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
    chi2 = float(np.sum((mean-M)**2/var))
    return M, sig, chi2, len(ti)-1

print("L1 F^2 0++ (Nf2 g8, t0=1) -- nops2 (rebase subspace) scan; 0++ = state nops2-2, fit [0.2,0.6]")
plt.figure(figsize=(7.4, 5.0))
ax = plt.gca()
ax.axhline(2*np.sqrt(2), ls="--", color="0.35", lw=1.1)
ax.text(0.22, 2*np.sqrt(2)+0.03, r"two-photon $2\sqrt{2}$", ha="left", va="bottom", color="0.35", fontsize=9)
ax.axvspan(rng[0], rng[1], color="0.85", alpha=0.5, zorder=0)
for NO in (2, 3, 4):
    jk, at = load(f"{SC}/no_f2_L1_no{NO}.dat")
    s0 = NO - 2                      # lightest 0++
    col, mk = nstyle[NO]
    n = min(ntshow, jk.shape[1])
    tt = np.array([(i+1)*at for i in range(n)]) + ndodge[NO]
    mean, err = curve(jk, s0)
    ax.errorbar(tt, mean[:n], yerr=err[:n], color=col, marker=mk, ms=6, capsize=2, lw=1.3,
                ls="-", label=f"nops2={NO}  (0++ = state {s0})", zorder=3)
    M, e, c, d = fit(jk, s0, rng[0], rng[1], at)
    ax.axhspan(M-e, M+e, color=col, alpha=0.10, zorder=1)
    print(f"   nops2={NO}: M_0++ = {M:.4f}({e:.4f})  chi2/dof={c/max(d,1):.2f}")
    # heavier excited 0++ (state s0-1) for nops2>=3, thin dashed
    if NO >= 3:
        me, ee = curve(jk, s0-1)
        ax.plot(tt, me[:n], color=col, lw=0.9, ls="--", alpha=0.55, zorder=2)
ax.plot([], [], color="0.4", lw=0.9, ls="--", label="excited 0++ (heavier state)")
ax.set_xlim(0.05, ntshow*at+0.05)
ax.set_ylim(1.8, 5.6)
ax.set_xlabel(r"$t$  [$a_t$ units, $t = n\,a_t$]")
ax.set_ylabel(r"$\Delta_\mathrm{eff}(t)$  [1/$a_t$]")
ax.set_title(r"Effective mass, Nf2 $g^2$=8, L=1: $F^2$ 0++ -- rebase-subspace (nops2) overlay")
ax.legend(loc="upper right", frameon=True, fontsize=8.5)
ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("glue_F2_nops_effmass_L1_claude.png", dpi=140)
plt.close()
print("wrote glue_F2_nops_effmass_L1_claude.png")
