# Overlaps among $\{\sigma_{PS},\ \sigma_{FS},\ \sigma_{PS}^2,\ \sigma_{FS}^2\}$

**Date:** 2026-09-17. Session "Better understand the {1,1,1,1}". Companion to
`fs_furnishing_derivation_claude.md` (the GW collapse of the FS leg) and
`sigma2_single_meson_exclusion_claude.md` (the $\sigma^2\to$ single-meson exclusion).

This note records the full two-point overlap structure of the four scalar operators, with the
reasoning for every entry. The headline correction relative to an earlier mid-session claim: the
single-single cross $\langle\sigma_{PS}\,\sigma_{FS}\rangle$ **vanishes** (opposite parity), it is
**not** equal to the auto-correlator, and it is **not** $O(a)$ -- it is an exact cancellation at
every $t\neq0$.

## Operators

The scalar bilinears (qed3int Eq. 5.1-5.3), $\sigma = \eta^\dagger S\,\xi + \xi^\dagger \tilde S\,\eta$
with $S=1$ throughout:

$$
\sigma_{PS} = \eta^\dagger \xi + \xi^\dagger \eta ,
\qquad
\sigma_{FS} = \eta^\dagger \xi + \xi^\dagger \tilde S\,\eta ,
\qquad
\tilde S = -(1 - D_\text{ov}^\dagger) .
$$

"PS" = **parity symmetric** furnishing ($\tilde S=+1$); "FS" = the furnished leg
($\tilde S=-(1-D_\text{ov}^\dagger)$). The squared operators $\sigma_{PS}^2$, $\sigma_{FS}^2$ are the
same bilinears at a common point (the $\sigma^2$ / two-bilinear operators).

## Contraction rules

The only nonzero elementary contractions are $\xi$ with $\eta^\dagger$ and $\eta$ with $\xi^\dagger$:

$$
\langle \xi\,\eta^\dagger\rangle = D_\text{ov}^{-1} \equiv \tau
\qquad\text{(forward inverse)} ,
$$
$$
\langle \eta\,\xi^\dagger\rangle = D_\text{ov}^{-\dagger} \equiv \tau^\dagger
\qquad\text{(adjoint inverse)} .
$$

$\langle\xi\,\xi^\dagger\rangle=\langle\eta\,\eta^\dagger\rangle=0$, so within a two-point only
like-with-like terms survive (no $\eta^\dagger\xi \times \xi^\dagger\eta$ cross).

## Two exact identities

**(1) Ginsparg-Wilson.** For the overlap operator,

$$
D_\text{ov}^{-1} + D_\text{ov}^{-\dagger} = 1
\qquad\Longrightarrow\qquad
\tau^\dagger = 1 - \tau .
$$

The "$1$" is the identity $\delta_{xy}\,\delta_{tt'}$. For unequal times $t\neq0$ it drops, so

$$
\boxed{\ \tau^\dagger(t,0) = -\,\tau(t,0)\quad (t\neq0)\ }
$$

This is exact (not $O(a)$): it is how continuum parity is realized on the lattice in this sector.

**(2) FS leg collapse.** The furnished leg contracts $\eta$ against $\xi^\dagger$ = the adjoint
inverse $\tau^\dagger$, dressed by $\tilde S$; by GW it collapses to the forward inverse:

$$
\tilde S\,\tau^\dagger = -(1 - D_\text{ov}^\dagger)\,D_\text{ov}^{-\dagger}
= -\big(D_\text{ov}^{-\dagger} - 1\big)
= 1 - D_\text{ov}^{-\dagger}
= D_\text{ov}^{-1} = \tau .
$$

So wherever the FS furnishing $\tilde S$ sits on a leg, that leg's $\tau^\dagger$ becomes $\tau$.
(Full derivation in `fs_furnishing_derivation_claude.md`.)

## Parity

In the continuum the two bilinears have **opposite** parity: one of $\{\sigma_{PS},\sigma_{FS}\}$ is
parity-even and the other parity-odd. The squared operators are therefore both parity-**even**
(even$^2$ = even, odd$^2$ = even). This is the organizing principle for the table.

## The overlap table (connected two-points, $t\neq0$)

|          | $\sigma_{PS}$ | $\sigma_{FS}$ | $\sigma_{PS}^2$ | $\sigma_{FS}^2$ |
|----------|:---:|:---:|:---:|:---:|
| $\sigma_{PS}$   | $M$ | $0$ | $0$ | $0$ |
| $\sigma_{FS}$   | $0$ | $M$ | $0$ | $0$ |
| $\sigma_{PS}^2$ | $0$ | $0$ | $T$ | $T$ |
| $\sigma_{FS}^2$ | $0$ | $0$ | $T$ | $T$ |

$M$ = single-meson two-point (nonzero); $T$ = two-meson four-point (nonzero); $0$ = exact zero for
$t\neq0$. Parity block-diagonalizes into a single-bilinear block and a $\sigma^2$ block; within the
single block PS and FS are orthogonal; within the $\sigma^2$ block PS$^2$ and FS$^2$ mix.

## Derivation of each block

### Single-single: PS-PS, FS-FS, PS-FS

Write $\sigma_{PS}(t)=\eta^\dagger_t\xi_t + \xi^\dagger_t\eta_t$ and
$\sigma_{FS}(0)=\eta^\dagger_0\xi_0 + \xi^\dagger_0\tilde S\,\eta_0$. Keeping the surviving Wick pairings:

$$
\langle\sigma_{PS}\,\sigma_{PS}\rangle
= -\mathrm{Tr}\!\big[\tau(t,0)\,\tau(0,t)\big]
  -\mathrm{Tr}\!\big[\tau^\dagger(t,0)\,\tau^\dagger(0,t)\big] ,
$$
$$
\langle\sigma_{FS}\,\sigma_{FS}\rangle
= -\mathrm{Tr}\!\big[\tau(t,0)\,\tau(0,t)\big]
  -\mathrm{Tr}\!\big[\tau(t,0)\,\tau(0,t)\big] ,
$$
$$
\langle\sigma_{PS}\,\sigma_{FS}\rangle
= -\mathrm{Tr}\!\big[\tau(t,0)\,\tau(0,t)\big]
  -\mathrm{Tr}\!\big[\tau^\dagger(t,0)\,\tau(0,t)\big] .
$$

(In FS-FS both furnished legs collapse $\tau^\dagger\to\tau$; in PS-FS only the one furnished leg
collapses.) Now apply $\tau^\dagger(t,0)=-\tau(t,0)$ for $t\neq0$:

$$
\langle\sigma_{PS}\,\sigma_{PS}\rangle
= -2\,\mathrm{Tr}[\tau(t,0)\tau(0,t)] \equiv M \neq 0 ,
$$
$$
\langle\sigma_{FS}\,\sigma_{FS}\rangle
= -2\,\mathrm{Tr}[\tau(t,0)\tau(0,t)] = M
\qquad(\text{same magnitude as PS-PS, hence bit-identical}) ,
$$
$$
\langle\sigma_{PS}\,\sigma_{FS}\rangle
= -\mathrm{Tr}[\tau\tau] + \mathrm{Tr}[\tau\tau] = 0 .
$$

**Reading.** PS and FS create the same scalar meson with the same $|$coupling$|$ (equal
auto-correlators), but they are **orthogonal** -- opposite parity. The PS-FS cancellation is exact at
every $t\neq0$ (it rides on the exact identity $\tau^\dagger=-\tau$ off-diagonal), so it is not an
$O(a)$ effect. $\sigma_{FS}$ is **not** the same operator as $\sigma_{PS}$; only their two-point
magnitudes coincide.

### Single-double (two-to-one): all zero

The four cross-grade overlaps $\langle\sigma_X\,\sigma_Y^2\rangle$ (one bilinear at the sink, two at
the source = a three-propagator triangle) all vanish, by two complementary mechanisms:

- **Parity.** $\langle\sigma_{FS}\,\sigma_{PS}^2\rangle$ and $\langle\sigma_{FS}\,\sigma_{FS}^2\rangle$
  are odd $\times$ even $\to 0$.
- **GW anti-hermiticity** (corrected — *not* $\sigma_3$-hermiticity; 2+1D 2-component has no chirality).
  The parity-allowed ones, $\langle\sigma_{PS}\,\sigma_{PS}^2\rangle$ and $\langle\sigma_{PS}\,\sigma_{FS}^2\rangle$
  (even $\times$ even), vanish because the normal-ordered propagator $M=D_\text{ov}^{-1}-\tfrac12$ is
  **anti-hermitian by GW** ($M^\dagger=-M$): a closed loop of $n$ $M$'s with hermitian vertices obeys
  $T^*=(-1)^n T$, so the **3-propagator** triangle (odd) is purely imaginary $\Rightarrow \mathrm{Re}\,T=0$
  (the physical correlator vanishes). This holds **config-by-config in the interacting theory**; the free case
  additionally kills the imaginary part site-by-site (the free $\sigma_3$-hermiticity). This is the
  $\sigma^2\to$ single-meson exclusion; full derivation in `gw_antiherm_exclusion_mechanism_claude.md`,
  summary in `sigma2_single_meson_exclusion_claude.md`.

So there is **no two-to-one overlap** in any PS/FS combination.

### Double-double: PS$^2$-PS$^2$, FS$^2$-FS$^2$, PS$^2$-FS$^2$

Both $\sigma^2$ operators are parity-even, so nothing in this block is parity-forbidden. The
auto-correlators are the two-meson four-point $T\neq0$, and $\langle\sigma_{FS}^2\,\sigma_{FS}^2\rangle
= \langle\sigma_{PS}^2\,\sigma_{PS}^2\rangle$ (the PP $==$ FF four-point equality, established
bit-identical, following from the leg-level $\tau^\dagger=-\tau$ collapse). The off-diagonal
$\langle\sigma_{PS}^2\,\sigma_{FS}^2\rangle$ is **nonzero** (even $\times$ even, not protected by any
trace identity here) -- the one surviving off-diagonal in the whole table.

### Why the four-point reduces to a single leg (audit)

Two facts make the four-point additive in a single leg object:

1. **S and S̃ never mix on a loop** (qed3int Eq. 5.5): a given fermion loop carries either all $S$-part
   legs or all $\tilde S$-part legs, so $\langle\sigma^2\sigma^2\rangle = \langle S\text{-part}\rangle +
   \langle\tilde S\text{-part}\rangle$ with a *single* leg type per part (no mixed-leg diagrams).
2. **Every assembled four-point diagram is homogeneous degree 4 in the leg** (a bilinear four-point has
   exactly four propagator lines), so it is invariant under a global $\text{leg}\to-\text{leg}$:
   $(-1)^4=+1$.

Combined with $\tau^\dagger=\delta-\tau$ (which gives $\tau^\dagger=-\tau$ off-diagonal *and*
$\tau^\dagger-\tfrac12=-(\tau-\tfrac12)$ at equal time), fact 2 gives $G_{10}[\tau^\dagger]=G_{10}[\tau]$.
Hence:

$$
\text{PS.PS} = G_{10}[\tau] + G_{10}[\tau^\dagger] = 2\,G_{10}[\tau],
\qquad
\text{FS.FS} = G_{10}[\tau] + G_{10}[\tilde S\,\tau^\dagger] = G_{10}[\tau] + G_{10}[\tau] = 2\,G_{10}[\tau].
$$

**Numerical audit** (free L1, `four_point_audit`, run 2026-09-17): with the honest adjoint
$\tau^\dagger=\delta-\tau$ (not the $-\tfrac12$ contact shortcut),

- each of the 10 assembled diagrams $G_{10,i}[\tau^\dagger]/G_{10,i}[\tau] = 1.000000$ (all even-degree),
- $G_{10}[\tau^\dagger] == G_{10}[\tau]$ to $3\times10^{-10}$,
- honest PS.PS $(G_{10}[\tau]+G_{10}[\tau^\dagger]) == 2\,G_{10}[\tau]$ to $2\times10^{-10}$,
- the artifact $G_{10}[-\tau_\text{gw}] = 0.041\times G_{10}[\tau]$ (the $\tau_\text{gw}$ bug).

Contrast with the PS-FS **2-point** cross, whose surviving diagram is degree 1 in $\tau^\dagger$ (odd) --
there the sign does *not* square away and the two terms cancel to $0$. Even-degree (four-point auto) $\to$
equal; odd-relative-degree (two-point cross) $\to$ cancels. Both are exact consequences of the same
$\tau^\dagger=-\tau$ identity.

## The $\tau_\text{gw}$ artifact (why an earlier note said PS-FS $\neq 0$)

Using `peram/tau_gw` $= V^\dagger(1-D_\text{ov}^\dagger)D_\text{ov}^{-1}V$ as the FS leg (the "$-\tau'$"
leg) is **wrong**: it carries the bare $D_\text{ov}^\dagger$ through the *forward* inverse and does not
satisfy $\tau^\dagger=-\tau$ off-diagonal. It therefore spoils the exact PS-FS cancellation and fakes a
nonzero cross (an earlier measurement read $R_{PF}\to -1$). The correct FS leg is the adjoint inverse
dressed by $\tilde S$, which collapses to plain $\tau$; with it, $\langle\sigma_{PS}\,\sigma_{FS}\rangle=0$
as required by parity. See the FIX SWEEP in `sigma2_single_meson_exclusion_claude.md`: FS leg
$= \tau$, never $-\tau_\text{gw}$.

## Summary

- Parity block-diagonalizes $\{\sigma_{PS},\sigma_{FS},\sigma_{PS}^2,\sigma_{FS}^2\}$ into a
  single-bilinear block and a $\sigma^2$ block.
- $\langle\sigma_{PS}\,\sigma_{PS}\rangle = \langle\sigma_{FS}\,\sigma_{FS}\rangle = M$ (same meson,
  equal magnitude); $\langle\sigma_{PS}\,\sigma_{FS}\rangle = 0$ **exactly** for $t\neq0$ (opposite
  parity), not $O(a)$.
- All four single-double (two-to-one) overlaps vanish -- by parity or by
  $\mathrm{Tr}[\Gamma_\text{even}GGG]=0$.
- The only off-diagonal survivor is $\langle\sigma_{PS}^2\,\sigma_{FS}^2\rangle\neq0$ (both
  parity-even).
