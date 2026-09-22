# Connected two-scalar $\sigma\sigma$ level -- GEVP benchmark

**Benchmark (what the $\{F^2,\sigma\sigma\}$ GEVP must reproduce in the $\sigma\sigma$ channel):**
$$
\boxed{\,a_t\,m_{\sigma\sigma}^\text{conn} \approx 1.0\quad\text{(dimensionless lattice units)}\,}
$$
Nf2 $g^2{=}0.5$ L1, $a_t{=}0.2$, massless, 400 cfg, distillation $N_v{=}24$ (exact). Code `two_meson_connected_claude.py`.

## Definition (fully-connected four-point cumulant)

$\sigma = \mathrm{PS}-\tfrac12$ (local scalar density, GW contact removed). Four points $x,y$ at time $t$, $z,w$ at $0$:
$$
\langle\sigma\sigma\sigma\sigma\rangle_c = \langle O(t)O(0)\rangle - \langle\sigma\sigma\rangle\langle\sigma\sigma\rangle - (\text{2 perms})
= \langle O(t)O(0)\rangle - \tfrac12\langle O\rangle^2 - 2\,\langle C_S(t)\rangle^2 ,
$$
with $O=\sigma^2$ (the two-meson operator), $C_S(t)=\langle\sigma(t)\sigma(0)\rangle$ the single-$\sigma$ two-point.
- **Vacuum term $\tfrac12\langle O\rangle^2$** (the $(xy)(zw)$ pairing): coefficient $\tfrac12$ from the $S_4{+}\tilde S_4$ flavor structure. Verified: $\mathrm{FP}(dt{\sim}30)=1.937{\times}10^{-4}$ vs $\tfrac12\langle O\rangle^2=1.940{\times}10^{-4}$ (0.2%).
- **Free two-$\sigma$ term $2\langle C_S(t)\rangle^2$** (the $(xz)(yw)+(xw)(yz)$ pairings): removes the two free single-$\sigma$ propagators (the scattering piece). Coefficient carries a flavor-factor caveat, but the vacuum-only (V0) and full (V1) effmasses agree to $\lesssim0.05$.

## Result and interpretation

- Connected effmass $a_t m \approx 0.9$–$1.0$ over the usable window dt$\,[4,7]$; signal dies into noise by dt$\sim$8 (no clean plateau -- still some excited-state weight, so $\sim1.0$ is approximate / upper-ish).
- This sits **above** the free two-$\sigma$ threshold $2m_\sigma\approx0.67$ and **above** the 0++ glueball $a_t m_{F^2}=0.616$. Removing the free piece *raises* the effmass ($0.92\to0.99$), i.e. slightly repulsive -- **no bound state below threshold**.
- Physics: the $\sigma\sigma$ channel is a **two-$\sigma$ scattering state** (discrete two-particle level on $S^2\times\mathbb{R}$), not a distinct light resonance.

## Reference scales (same ensemble, $a_t m$ lattice units)

| quantity | $a_t m$ | $\Delta/\Delta_A$ |
|---|---|---|
| axial $\ell{=}1$ (yardstick $\Delta_A$) | $0.357$ | $1$ |
| single $\sigma$ ($C_S$ plateau) | $\sim0.35$ | $\sim0.98$ |
| free two-$\sigma$ $2m_\sigma$ | $\sim0.67$ | $\sim1.9$ |
| 0++ glueball $F^2$ | $0.616$ | $1.72$ |
| **connected $\sigma\sigma$ (this benchmark)** | $\sim1.0$ | $\sim2.8$ |

Note $\Delta/\Delta_A$ uses physical $m=a_t m/a_t$ over the axial's physical $m$ (equivalently the $a_t m$ ratio); $a_t m_A=0.357$ is the redo record's $1.787$ times $a_t$ (their records are physical $m=\mathrm{acosh}/a_t$; ours are lattice $a_t m$).

## Next

Coupled $\{F^2,\sigma\sigma\}$ GEVP: the $\sigma\sigma$ operator should land near this $a_t m\approx1.0$; the question is whether the 0++ ($0.616$) mixes in as a distinct level. Also: derive the free-two-$\sigma$ flavor coefficient rigorously to sharpen V1.
