# FS furnishing: correct contraction is plain $\tau$ (GW collapse) -- $\tau_{gw}$ is the WRONG leg

## The question
Does $\sigma_{FS}^2$ couple to a single meson? Earlier I measured O(a) (PS$^2\to$FS) and O(1) (FS$^2\to$PS/FS)
using the stored `peram/tau_gw` as the FS leg. NM's challenge: the two-flavor propagator adjoint
($\langle\eta\,\xi^H\rangle=D_{\rm ov}^{-\dagger}$, NOT $D_{\rm ov}^{-1}$) means those couplings may be a
$D^{-1}$-vs-$D^{-\dagger}$ build artifact. This note derives the correct FS leg from the definition.

## Definition (qed3int Eq 5.1-5.3)
$$\sigma_{PS/FS} = \eta^\dagger S\,\xi + \xi^\dagger \tilde S\,\eta,\qquad S=1,\quad
\tilde S = 1\ (\text{PS})\ /\ -(1-D_{\rm ov}^\dagger)\ (\text{FS}).$$
$\tilde S$ is a **vertex operator**, not a propagator. $S,\tilde S$ never mix on a loop (Eq 5.5):
$\langle\sigma\sigma\rangle=\langle S\text{-part}\rangle+\langle\tilde S\text{-part}\rangle$.

## Propagators (two-flavor / four-component) + GW (IV.17)
$$\langle\xi\,\eta^H\rangle = D_{\rm ov}^{-1}\ (=\tau),\qquad
  \langle\eta\,\xi^H\rangle = D_{\rm ov}^{-\dagger}\ (=D_{\rm ov}^{-H}),\qquad
  D_{\rm ov}^{-1}+D_{\rm ov}^{-\dagger}=1 .$$
So the ADJOINT inverse is bounded: $D_{\rm ov}^{-\dagger}=1-D_{\rm ov}^{-1}=\delta-\tau$ (in the distillation
basis, complete: $=I-\tau$), available for FREE from `peram/tau`, NO extra solve.

## The collapse (the whole point)
- **$S$-part** (both vertices $S=1$, fields $\eta^\dagger\xi$): legs contract $\xi\!\leftrightarrow\!\eta^H$
  $=D_{\rm ov}^{-1}=\tau$. Plain. Identical to PS.
- **$\tilde S$-part** (both vertices $\tilde S=-(1-D_{\rm ov}^\dagger)$, fields $\xi^\dagger\eta$): each leg is the
  vertex operator $\tilde S$ acting through the ADJOINT propagator $\langle\eta\,\xi^H\rangle=D_{\rm ov}^{-\dagger}$:
$$\tilde S\cdot D_{\rm ov}^{-\dagger} = -(1-D_{\rm ov}^\dagger)\,D_{\rm ov}^{-\dagger}
  = -\big(D_{\rm ov}^{-\dagger} - D_{\rm ov}^\dagger D_{\rm ov}^{-\dagger}\big)
  = -\big(D_{\rm ov}^{-\dagger} - 1\big) = 1 - D_{\rm ov}^{-\dagger} = D_{\rm ov}^{-1} = \tau .$$
  (used $D_{\rm ov}^\dagger D_{\rm ov}^{-\dagger}=1$ and GW $1-D_{\rm ov}^{-\dagger}=D_{\rm ov}^{-1}$.)

**So the FS leg collapses to plain $\tau$.** $\tilde S$-part $=$ $S$-part, hence $\sigma_{FS}$ computes
**identically to $\sigma_{PS}$**: $\langle\sigma_{FS}\sigma_{FS}\rangle=\langle\sigma_{PS}\sigma_{PS}\rangle$
(this is the known GW "PS$=$FS", and matches the bit-identical PP$==$FF cache). Both are the clean
$(\bar\psi\psi)^2$, **protected**, and neither couples to a single meson.

## Why my earlier result was an artifact
The stored `peram/tau_gw` $= V^\dagger(1-D_{\rm ov}^\dagger)D_{\rm ov}^{-1}V$ uses the **forward** $D_{\rm ov}^{-1}$,
not the adjoint $D_{\rm ov}^{-\dagger}$. By GW $(1-D_{\rm ov}^\dagger)D_{\rm ov}^{-1}=1+D_{\rm ov}^{-1}-D_{\rm ov}^\dagger$
-- it carries the **bare $D_{\rm ov}^\dagger$** operator (O(1), unbounded). Using $-\tau_{gw}$ as the FS leg is
therefore NOT the Eq-5.1-5.3 contraction; the bare-$D_{\rm ov}^\dagger$ piece manufactured the spurious O(a)
(PS$^2\to$FS) and O(1) (FS$^2\to$PS/FS) single-meson couplings. **They are not physical.**

## Which codebase object is which (Fin, confirmed independently)
- **flavfac $FF$** (`sigma2_flavorgeom` cache: AblkS$=\tau$ $+$ per-loop $(1+(-1)^{n_{FS}})$): this **IS the correct
  $\sigma_{FS}^2$** -- the per-loop sign IS the collapse-to-$\tau$, so $==$ PS $==$ PP. Not a "lesser protected
  version." Fin's $\sigma^2$-$F^2$ 6x6 GEVP uses this cache -> UNAFFECTED.
- **production `four_point` FS.FS** $=S_4[\tau]+S_4[-\tau_{gw}]$ (`distill_contract:284`): carries the artifact
  (measured $\approx0.52\,$PS.PS, flat; should $==$PS.PS $=2S_4[\tau]$). This is the buggy object.
- **PRODUCTION-SCOPE FLAG:** anything using `tau_gw`/$(1-D_{\rm ov}^\dagger)$ at higher-than-2pt order likely
  carries the same bare-$D_{\rm ov}^\dagger$. The FS 2pt survived ($\langle\sigma_{FS}\sigma_{FS}\rangle\to0.378$)
  because the 2pt is insensitive; $\sigma^2$/triangles are where it bites. (Fin is raising a `tau_gw`-usage audit.)

## Correct implementation (NM: use $\tau$ only)
The FS leg is just $\tau$ (the collapse). Equivalently, if one keeps the $\tilde S$ vertex explicit, the leg is
the adjoint propagator $D_{\rm ov}^{-\dagger}=\delta-\tau$ dressed by $\tilde S$ -- but that product IS $\tau$, so
**use `peram/tau` and never `peram/tau_gw`** for these objects. `tau_gw` (with its bare $D_{\rm ov}^\dagger$) is
only correct where the literal operator $(1-D_{\rm ov}^\dagger)$ acting on the forward inverse is what is wanted.

## Consequence (corrects `sigma2_single_meson_exclusion_claude.md`)
- $\sigma_{PS}^2$ and $\sigma_{FS}^2$ are the SAME correlator (GW), both clean $(\bar\psi\psi)^2$, **protected**.
- $\sigma^2$ (PS or FS) does **NOT** couple to any single meson. The earlier O(a)/O(1) FS couplings were
  the `tau_gw` ($D^{-1}$-vs-$D^{-\dagger}$) artifact, NOT physics.
- The only exact statement that stands: $\langle\sigma_{00}\sigma^2\rangle=0$ (GW anti-hermiticity of $M=D_{\rm ov}^{-1}-\tfrac12$, odd 3-loop $\to$ Re$=0$; free-only $\sigma_3$-herm; NOT interacting $\sigma_3$-herm -- see `gw_antiherm_exclusion_mechanism_claude.md`), and
  PS$==$FS. There are NOT "two inequivalent constructions" -- there is one correct contraction ($\tau$), and a
  wrong leg ($\tau_{gw}$).

Refs: qed3int Eq 5.1-5.5 (furnishing), IV.17/1.36 (GW $D_{\rm ov}^{-1}+D_{\rm ov}^{-\dagger}=1$);
`distill_peram_claude.cu` (tau/tau_gw build); T2b validation $\bar\tau=\delta-\tau=D_{\rm ov}^{-\dagger}$.
