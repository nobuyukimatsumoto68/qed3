#!/usr/bin/env python3
# cheby_filter_plot_claude.py
# _claude: plot the EXACT Chebyshev filter used by the IRL run, p(mu) = T_n(q(mu)), applied to
# A = Dov^dag Dov (eigenvalue mu in [0,4]).  Header convention (lanczos_claude.h):
#   q(mu) = (2 mu - (alpha^2 + beta^2)) / (alpha^2 - beta^2)     [alpha,beta are SINGULAR-value edges]
#   T_n via the three-term recurrence T_0=1, T_1=q, T_{k+1}=2 q T_k - T_{k-1}  (exactly the header's loop)
# Shows whether the wanted low band (mu < beta^2) is amplified relative to the bulk.

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

alpha = 2.0        # header arg; alpha^2 = 4 = max eigenvalue of Dov^2 (spectral bound)
beta = 0.5         # header arg; beta^2 = 0.25 = wanted-band edge (mu < 0.25, |z| < 0.5 = near-zero tail)
degs = [12, 8, 6]  # compare Chebyshev orders (gain at mu->0 grows as cosh(n*arccosh|q(0)|))

a2 = alpha * alpha
b2 = beta * beta

mu = np.linspace(0.0, 4.0, 4000)
q = (2.0 * mu - (a2 + b2)) / (a2 - b2)


def cheb(qq, n):
    # T_n(qq) via the same three-term recurrence the header runs (handles |q|>1 -> grows).
    Tprev = np.ones_like(qq)     # T_0
    Tcur = qq * np.ones_like(qq) # T_1
    for k in range(2, n + 1):
        Tprev, Tcur = Tcur, 2.0 * qq * Tcur - Tprev
    return Tcur


fig, ax = plt.subplots(figsize=(8.5, 5.5))
styles = [("C0", "-"), ("C1", "--"), ("C2", "-.")]
mu_expected = 0.05
q0 = (2.0 * mu_expected - (a2 + b2)) / (a2 - b2)
for (n, (col, ls)) in zip(degs, styles):
    absp = np.abs(cheb(q, n))
    gain0 = float(np.abs(cheb(np.array([q0]), n))[0])
    ax.semilogy(mu, absp, color=col, ls=ls, lw=1.6,
                label=r"$n=%d$  (gain@$\mu{=}0.05$ = %.1e)" % (n, gain0))

ax.axvline(b2, color="C3", ls=":", lw=1.2, label=r"wanted edge $\beta^2=%.2f$" % b2)
ax.axvline(a2, color="0.4", ls=":", lw=1.2, label=r"spectral bound $\alpha^2=%.2f$" % a2)

ax.set_xlabel(r"$\mu$  (eigenvalue of $A=D_{ov}^\dagger D_{ov}$, range $[0,4]$)")
ax.set_ylabel(r"$|T_n(q(\mu))|$  (filter gain, log scale)")
ax.set_title(r"Chebyshev filter vs order:  $\alpha=%.1f,\ \beta=%.1f$" % (alpha, beta))
ax.set_xlim(0.0, 4.0)
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="upper center", fontsize=8, framealpha=0.9)

fig.tight_layout()
fig.savefig("cheby_filter_claude.png", dpi=130)
print("wrote cheby_filter_claude.png")
for n in degs:
    g = float(np.abs(cheb(np.array([q0]), n))[0])
    print("deg=%2d : gain(mu=0.05) = %.3e" % (n, g))
