#!/usr/bin/env python3
# A/B effective mass: spatial-only flow (per-timeslice, get_spatial) vs full-3D flow (get_force),
# Nf2 g8, same configs. 4 panels-as-figures (F/F^2 x L1/L2). F = unit m-average of the per-m ground;
# F^2 = state0 (0++). Bands = diagonal-chi2 const fit (F [0.2,0.8], F^2 [0.2,0.6]). Anchors sqrt2 / 2sqrt2.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SC = "/tmp/claude-1000/-mnt-barracuda22-qed3/786b05aa-b873-4bf5-9b42-33d676565057/scratchpad"
fstyle = {"spatial": ("#1f77b4", "o", "spatial-only flow (per-timeslice)"),
          "full": ("#d62728", "s", "full-3D flow (get_force)")}
fdodge = {"spatial": -0.012, "full": 0.012}

def load(fn):
    h = open(fn).readline().split()
    nb, nt, ns, at = int(h[1]), int(h[2]), int(h[3]), float(h[4])
    return np.loadtxt(fn).reshape(nb, nt, ns), at

def curve(jk, chan):
    nb = jk.shape[0]
    c = jk.mean(axis=2) if chan == "F" else jk[:, :, 0]
    mean = c.mean(0)
    err = np.sqrt((nb-1)/nb*np.sum((c-mean)**2, 0))
    return mean, err

def fitband(jk, chan, tlo, thi, at):
    nb, nt, ns = jk.shape
    ti = [t for t in range(nt) if tlo-1e-9 <= (t+1)*at <= thi+1e-9]
    if chan == "F":
        Mbs = []
        sigs = []
        for s in range(ns):
            d = jk[:, ti, s]
            var = (nb-1)/nb*np.sum((d-d.mean(0))**2, 0)
            w = 1.0/var
            w /= w.sum()
            Mbs.append(d @ w)
            sigs.append(np.sqrt((nb-1)/nb*np.sum((Mbs[-1]-Mbs[-1].mean())**2)))
        Mbs = np.array(Mbs)
        u = np.array([1.0/s**2 for s in sigs])
        u /= u.sum()
        Mc = (u[:, None]*Mbs).sum(0)
        M = Mc.mean()
        e = np.sqrt((nb-1)/nb*np.sum((Mc-M)**2))
    else:
        d = jk[:, ti, 0]
        var = (nb-1)/nb*np.sum((d-d.mean(0))**2, 0)
        w = 1.0/var
        w /= w.sum()
        Mb = d @ w
        M = Mb.mean()
        e = np.sqrt((nb-1)/nb*np.sum((Mb-M)**2))
    mean = jk[:, ti, 0].mean(0) if chan != "F" else None
    r = (jk[:, :, :].mean(0)[ti].mean(1) if False else None)
    # chi2 (diag) of the fitted state
    s0 = 0
    dd = jk[:, ti, s0]
    mn = dd.mean(0)
    vv = (nb-1)/nb*np.sum((dd-mn)**2, 0)
    chi2 = float(np.sum((mn-M)**2/vv))/max(len(ti)-1, 1)
    return M, e, chi2

panels = [("F", 1, (0.2, 0.8), np.sqrt(2), r"$\sqrt{2}$", (1.0, 1.6), 8),
          ("F", 2, (0.2, 0.8), np.sqrt(2), r"$\sqrt{2}$", (0.8, 1.5), 8),
          ("F2", 1, (0.2, 0.6), 2*np.sqrt(2), r"$2\sqrt{2}$", (1.8, 3.4), 7),
          ("F2", 2, (0.2, 0.6), 2*np.sqrt(2), r"$2\sqrt{2}$", (1.5, 4.2), 7)]

for chan, L, rng, anchor, alab, ylim, ntshow in panels:
    plt.figure(figsize=(7.2, 5.0))
    ax = plt.gca()
    ax.axhline(anchor, ls="--", color="0.35", lw=1.1)
    ax.text(0.22, anchor, "free cont. "+alab, ha="left", va="bottom", color="0.35", fontsize=9)
    ax.axvspan(rng[0], rng[1], color="0.85", alpha=0.5, zorder=0)
    for flow in ("spatial", "full"):
        jk, at = load(f"{SC}/ab_{chan}_L{L}_{flow}.dat")
        col, mk, lab = fstyle[flow]
        n = min(ntshow, jk.shape[1])
        tt = np.array([(i+1)*at for i in range(n)]) + fdodge[flow]
        mean, err = curve(jk, chan)
        M, e, ch = fitband(jk, chan, rng[0], rng[1], at)
        ax.errorbar(tt, mean[:n], yerr=err[:n], color=col, marker=mk, ms=6, capsize=2, lw=1.5,
                    ls="-", label=f"{lab}:  M={M:.3f}({e:.3f}), $\\chi^2$/dof={ch:.2f}", zorder=3)
        ax.axhspan(M-e, M+e, color=col, alpha=0.12, zorder=1)
    ax.set_xlim(0.05, ntshow*at+0.05)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$t$  [$a_t$ units, $t = n\,a_t$]")
    ax.set_ylabel(r"$\Delta_\mathrm{eff}(t)$  [1/$a_t$]")
    ttl = "linear-$F$ $\\ell=1$" if chan == "F" else "$F^2$ 0++"
    ax.set_title(f"Flow A/B, Nf2 $g^2$=8, L={L}: {ttl}")
    ax.legend(loc="upper right", frameon=True, fontsize=8.5)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    out = f"glue_{chan}_flowAB_L{L}_claude.png"
    plt.savefig(out, dpi=140)
    plt.close()
    print("wrote", out)
