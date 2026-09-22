# Distillation / LapH for the two-meson ($\sigma\sigma$) four-point in QED3

Working reference for using distillation to compute the connected $\sigma\sigma$ four-point (the ten
diagrams A--J of `qed3int_v3-4.pdf` Ch.5), which is otherwise blocked by the all-to-all propagator
problem in $S_S/T_S$ (the internal $t\to t$ return, diagrams A, B). Also the noise-mitigation tool for
the disconnected loops. Specialized to our setup: $S^2\times\mathbb{R}$, $U(1)$ gauge, MASSLESS overlap
$D_{ov}$, 2-component GW fermions, NO $\gamma_5$-hermiticity.

**References (mandatory citation):** M. Peardon et al. (Hadron Spectrum), *Distillation*, arXiv:0905.2160;
C. Morningstar et al., *Stochastic LapH*, arXiv:1104.3870; all-to-all/dilution J. Foley et al.,
hep-lat/0505023; low-mode averaging Giusti hep-lat/0402002, DeGrand-Schaefer hep-lat/0401011.

---

## 1. Why distillation for our case

The $\sigma\sigma$ four-point needs the fermion propagator between ALL timeslice pairs (the internal
$\Pi_t$ return in $S_S/T_S$), which point/wall sources give only at $O(N_t)$ solves. Distillation solves
BOTH problems at once:

1. **All-to-all.** The perambulator (below) is the propagator between every timeslice pair, computed once
   per config. Every diagram -- including the internal returns -- becomes a trace over perambulators and
   vertex matrices, with NO per-$t$ and NO per-diagram inversion.
2. **Noise.** Full distillation is EXACT in the smeared subspace (no stochastic noise), so the notoriously
   noisy disconnected loops $D_S,D'_S$ and the four-point carry only GAUGE variance. And the smeared
   operators overlap the low-lying $0^{++}$ states much better -> earlier, flatter plateaus.

Smearing does NOT bias the mixing conclusion: it changes operator overlaps, not the states. GEVP
eigenvalues (masses) and whether the light branch carries both $F^2$ and $\sigma\sigma$ content are
physical and smearing-independent.

---

## 2. The distillation subspace on $S^2\times\mathbb{R}$

Spatial sites at timeslice $t$: $x=1,\dots,N_s$, $N_s=10L^2+2$ (the fixed $S^2$ triangulation). Spin
$\alpha=1,2$ (2-component GW; NO color in $U(1)$).

Build the **gauge-covariant spatial graph Laplacian** at each timeslice $t$ from the $U(1)$ spatial links
$U_{xy}(t)$ on the $S^2$ mesh:
$$
(-\Delta_t)_{xy}=\kappa_x\,\delta_{xy}-\!\!\sum_{\langle xy\rangle}\!c_{xy}\,U_{xy}(t),
$$
($c_{xy}$ = geometric edge weights of the simplicial mesh, $\kappa_x=\sum_y c_{xy}$; acts on the SITE
index only, identity on spin). Take its lowest $N_v$ eigenvectors
$$
V(t)=[\,v_1(t),\dots,v_{N_v}(t)\,],\qquad v_k(t)\in\mathbb{C}^{N_s},\qquad -\Delta_t\,v_k=\lambda_k v_k,\ \lambda_1\le\dots\le\lambda_{N_v}.
$$
The **smearing (distillation) operator** at $t$ is the projector onto that subspace, identity on spin:
$$
\Box(t)=V(t)V^\dagger(t)\quad(N_s\times N_s,\ \text{per spin}).
$$
$N_v$ sets the smearing radius; on the fixed $S^2$ mesh $N_v$ is modest (tune by the effective smearing
profile / plateau quality). $V(t)$ is per-config, per-timeslice (gauge-covariant).

---

## 3. Perambulators (overlap, NO $\gamma_5$-herm)

The **perambulator** is the all-to-all propagator projected into the distillation subspace, carrying the
$2\times2$ spin structure of $D_{ov}$:
$$
\tau^{\alpha\beta}_{kl}(t',t)\ \equiv\ v_k^\dagger(t')\,\big[D_{ov}^{-1}\big]^{\alpha\beta}(t',t)\,v_l(t).
$$
Compute it by solving, for each source $(l,\beta,t)$, the smeared source $S_{l,\beta}(t)=v_l(t)\otimes e_\beta$
(supported on timeslice $t$):
$$
D_{ov}\,\psi=S_{l,\beta}(t)\ \Rightarrow\ \psi=D_{ov}^{-1}S_{l,\beta}(t),\qquad
\tau^{\alpha\beta}_{kl}(t',t)=v_k^\dagger(t')\,\psi^\alpha(t').
$$
Cost: $N_v\times 2\times N_t$ forward solves (source at every timeslice for the all-to-all).

**No $\gamma_5$-herm -> also need the BACKWARD perambulator.** The meson legs use both $D_{ov}^{-1}$ and
$D_{ov}^{-\dagger}$ (there is no $\gamma_5$-herm to relate them; verified via the normal-equation route
$D_{ov}^{-\dagger}b=D_{ov}(D_{ov}^\dagger D_{ov})^{-1}b$, template `jj_local_ylm_scalar_conn_stoch_claude.cu`):
$$
\bar\tau^{\alpha\beta}_{kl}(t',t)\equiv v_k^\dagger(t')\,\big[D_{ov}^{-\dagger}\big]^{\alpha\beta}(t',t)\,v_l(t),
\qquad\text{another }N_v\times2\times N_t\text{ solves (same cost as }\tau).
$$
Both come from the SAME $V(t)$ and the SAME $D_{ov}^\dagger D_{ov}$ CG (differ only by which apply follows),
so fwd+bwd is 2x the solves, not more machinery.

---

## 4. Elementals (the $\ell{=}0$ scalar vertex)

The scalar density $\sigma=\eta^\dagger S\xi+\xi^\dagger\tilde S\eta$ ($S=1$; $\tilde S=1$ PS / $-(1-D_{ov}^\dagger)$
FS) becomes, projected into the subspace at timeslice $t$, an $N_v\times N_v$ (with spin) **elemental**:
$$
\Phi^{\alpha\beta}_{kl}(t)=v_k^\dagger(t)\,\big[W_0\,\Gamma\big]^{\alpha\beta}\,v_l(t),
$$
with $W_0=$ the $\ell{=}0$ area weight $A_x$ (spatial-sum, zero-momentum; higher $\ell$ carry $Y_{\ell m}$).
- **PS** ($S=1$): $\Gamma=\mathbb{1}$ -> $\Phi^{PS}(t)=V^\dagger(t)\,W_0\,V(t)$ (a smeared identity-density, cheap, no solve).
- **FS** ($\tilde S=-(1-D_{ov}^\dagger)$): the $(1-D_{ov}^\dagger)$ is a full operator, so it does NOT close
  inside one timeslice's $V(t)$. Fold it into a **generalized perambulator**: define
  $\tau'^{\alpha\beta}_{kl}(t',t)=v_k^\dagger(t')\,[(1-D_{ov}^\dagger)D_{ov}^{-1}]^{\alpha\beta}(t',t)\,v_l(t)$
  (one extra $D_{ov}^\dagger$ apply on $\psi$ before the $V^\dagger$ contraction -- no extra solve). FS lines
  use $\tau'$ where PS lines use $\tau$. (This mirrors the loop driver's $J_1$ vs $J_{1mD}$ split.)

---

## 5. The ten diagrams as perambulator traces

With $\{\tau,\bar\tau,\tau',\Phi\}$ (all $N_v\!\times\!N_v$ with spin, per config), every diagram is a
trace -- the internal $\Pi_t$ return is just the equal-time perambulator $\tau(t,t)$, an index contraction.
Schematically (indices/traces over $\{k,\alpha\}$; $\Phi$ at the vertices, $\tau/\bar\tau$ on the lines):

| factor | distillation trace | note |
|---|---|---|
| $D_S(t)$ (5.15) | $\mathrm{Tr}[\Phi(t)\,\tau(t,t)]$ | disc loop = equal-time perambulator |
| $D'_S(t)$ (5.16) | $\mathrm{Tr}[\Phi(t)\,\tau(t,t)\,\Phi(t)\,\tau(t,t)]$ | internal $\Pi_t\to\tau(t,t)$ |
| $C_S(0{\to}t)$ (5.17) | $\mathrm{Tr}[\Phi(0)\,\tau(0,t)\,\Phi(t)\,\bar\tau(t,0)]$ | meson corr, fwd + bwd legs |
| $V_S$ (5.18) | $\mathrm{Tr}[\Phi(0)\,\tau(0,0)\,\Phi(0)\,\tau(0,t)\,\Phi(t)\,\bar\tau(t,0)]$ | one $\Pi_t$ |
| $S_S$ (5.19) A | $\mathrm{Tr}[\Phi(0)\,\tau(0,0)\,\Phi(0)\,\tau(0,t)\,\Phi(t)\,\tau(t,t)\,\Phi(t)\,\bar\tau(t,0)]$ | TWO $\Pi_t\to\tau(t,t)$ |
| $T_S$ (5.20) B | $\mathrm{Tr}[\Phi(0)\,\tau(0,t)\,\Phi(t)\,\tau(t,0)\,\Phi(0)\,\tau(0,t)\,\Phi(t)\,\bar\tau(t,0)]$ | $(\Pi_t,\Pi_0,\Pi_t)$ order |

(Exact index/spin bookkeeping, the $\Pi_0$ vs $\Pi_t$ placement, and FS $\tau\to\tau'$ substitution per
vertex to be pinned at implementation against Ch.5; the point is every diagram is a finite product of the
precomputed $N_v\!\times\!N_v$ matrices -- no per-$t$ solve, connected & disconnected uniform.) The 10
diagrams then combine with the weights $\{4,2,4,4,2,1,4,1,1,1\}$; squares (E,H,I,J) still need the
distinct-source rule, but distillation being EXACT removes the stochastic self-variance entirely (the
$v/H$ bias of `two_meson_square_bias_note_claude.md` vanishes -- no hits).

---

## 6. The $F^2$ cross and the disc block

- **$F^2$ cross** $\langle F^2(t)\,\sigma\sigma(0)\rangle_c=\langle O_F(t)\,[D_S(0)^2+D'_S(0)]\rangle_c$:
  the scalar side is $D_S(0)$ and $D'_S(0)$ from the equal-time perambulator $\tau(0,0)$ (Sec. 5), times the
  gluonic $O_F(t)$ -- no new perambulator, no all-to-all. (This is the mixing-CHECK observable; already
  launched stochastically via the loop driver, but distillation would de-noise it.)
- **$\sigma\sigma$ disc self-correlator**: products of the single-loop 2pt $C_{xy}(t)=\langle J_x(t)J_y(0)\rangle$;
  each $J_x(t)=D_S^x(t)=\mathrm{Tr}[\Phi_x(t)\tau(t,t)]$ from the perambulator.

---

## 7. Cost and parameters

- **Eigensolve:** lowest $N_v$ eigenvectors of the $U(1)$-covariant $-\Delta_t$ on the $S^2$ mesh, per
  (config, timeslice). $N_s=10L^2+2$ small at L1/L2 -> cheap (dense or Lanczos; reuse `lanczos_claude.h`).
- **Perambulators:** $\tau$ and $\bar\tau$ each $N_v\times2\times N_t$ overlap solves per config
  ($=2\times$ that with the backward one). Generalized $\tau'$ (FS) = one extra $D_{ov}^\dagger$ apply on the
  same $\psi$, no extra solve.
- **Elementals:** $\Phi(t)=V^\dagger(t)W_0V(t)$ per timeslice, negligible.
- **Contractions:** each diagram = $O(N_t)$ products of $N_v\!\times\!N_v$ matrices, cheap offline.
- **Tunable:** $N_v$ (smearing radius / cost tradeoff); pick by plateau quality on a few configs. Start
  small ($N_v\sim$ tens) at L1.

---

## 8. Practical / implementation notes

- Overlap apply inside the perambulator solve = the Zolotarev sign function (inner multishift CG); this is
  the dominant cost, same as any overlap solve. Reuse the existing `overlap_wmass_claude.h` (mass 0) +
  `matpoly` machinery.
- Store $\tau,\bar\tau,\tau'$ per config (small: $N_v^2\,N_t^2\,\times$ spin complex) for offline diagram
  contraction and the coupled GEVP.
- Connection to LMA: full distillation here is EXACT in the subspace; if the low-Laplacian modes miss
  relevant high-mode content, the LMA completion (exact-low + stochastic-high remainder, Sec. 1) restores
  the LOCAL operator -- deferred, only if the smeared operator's overlap proves insufficient.
- This SUPERSEDES the two-one-end-trick sketch for A/B: distillation handles the internal $\Pi_t$ generally
  and de-noises, rather than a bespoke second-source trick per topology.

---

*Status: candidate method for the FULL $\sigma\sigma$ conn / coupled GEVP (deferred). The mixing CHECK
(loop driver, launched) does not need distillation. See `sigma_sigma_f2_mixing_impl_plan_claude.md`.*
