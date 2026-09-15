# 2PS $\to$ PS/FS three-point ($\langle\sigma^2\,\sigma\rangle$, the $2\to1$ triangle) -- plan

The one place FS is NOT redundant with PS: an ODD-vertex correlator. The cleanest is the $2\to1$ three-point
$\langle\sigma_{PS}^2(t)\,\sigma_X(s)\rangle$, $X\in\{PS,FS\}$ -- a single triangle (3-loop).

**All in the $\ell=0$ ($0^{++}$) channel.** Every operator here is Y00-projected / rotationally invariant:
$\sigma^2_{00}$ (verified couples only to $\ell=0$), $\sigma_{FS,00}$ (Y00), and the point operators $O_{2m}$,
$O_{1m}$ (rotational scalars via the $x$-sum with the rotation-covariant antipode/coincident pairing). So the
single-meson content below is the $\ell=0$ single meson ($m_{PS}$).

## Expected content = SINGLE-MESON tower (not the two-meson)
$$
\langle\sigma_{PS}^2(t)\,\sigma_X(s)\rangle = \sum_n \langle0|\sigma_{PS}^2|n\rangle\,\langle n|\sigma_X|0\rangle\,e^{-E_n(t-s)} .
$$
$\sigma_X$ is a BILINEAR (2-fermion). By quantum numbers nothing forbids $\langle{\rm two\text{-}meson}|\sigma_X|0\rangle$
(both are $0^{++}$; the second $q\bar q$ can come from the sea) -- so this is NOT a strict selection rule (an earlier
version wrongly said "unreachable by particle number"; corrected). What is true is that the CONNECTED (valence)
triangle contraction crosses the time cut with a single $q\bar q$, so it reaches the **single-meson tower**; the
two-meson enters only through DISCONNECTED / sea pieces, which are dynamically suppressed. Hence the (connected)
correlator is dominated by the $0^{++}$ single-meson tower: **ground $\to m_{PS}\approx0.32$** (and single-meson
excited states), NOT $2m_{PS}$. This is a DIAGNOSTIC -- seeing $m_{PS}$ confirms the cubic $2\to1$ coupling
$\langle0|\sigma^2|{\rm meson}\rangle$ (the fermion triangle) and the single-meson valence content.

## Diagram = the triangle (= diagram A / $C_{2A}$ structure), one FS leg for $X=FS$
3 $\sigma$ vertices (2 at $t$, 1 at $s$) -> 3 propagators -> one connected loop (triangle):
$$
T_X(t,s) = -\mathrm{Tr}\big[\,\Phi(t)\,\tilde\tau(t,t)\,\Phi(t)\,\text{leg}_X(t,s)\,\Phi_X(s)\,\text{leg}_X(s,t)\,\big] + (\text{sym}),
$$
one EQUAL-TIME leg $\tilde\tau(t,t)$ (improved, so any pure tadpole piece vanishes) + two off-diagonal legs.
- $X=PS$: all legs $\tau$ (forward improved). This is the all-PS triangle = the $C_{2A}$ / diagram-A object we
  already compute (the single-meson contamination of $\sigma^2$).
- $X=FS$: the $\sigma_{FS}$ vertex at $s$ is the $\tilde S$ one -> its two legs are the BACKWARD-improved
  ($-\tilde\tau$, from $(1-D_{ov}^\dagger)D^{-\dagger}=-D^{-1}$). Because the triangle is a lone 3-loop (odd),
  the $(-1)^3$ / furnishing SURVIVES (no tadpole to mask it): $T_{FS}\ne T_{PS}$. This is the nontrivial channel.

So $T_{PS}$ and $T_{FS}$ have the SAME spectrum (single-meson) but DIFFERENT overlaps -> different amplitudes,
same $m_{PS}$ plateau. Comparing them isolates the furnishing's effect on the $2\to1$ coupling.

## Computation
- Reuse the diagram-A / $C_{2A}$ machinery (`gevp_twosigma_interacting` has $C_{2A}=2\,$Tri; `diag_OA_sig2`).
  For $X=FS$ swap the $s$-vertex legs to backward-improved ($-\tilde\tau$ equal-time / $-\tau$ off-diagonal).
  Always the $1/2$-improved propagator. Translation-average $s$ over the window; the same triple subtraction
  (per-config t-sum + plateau) for any residual constant; config jackknife; log$|T|$ + effmass.
- Free-theory precedent: `two_meson_ps2_to_fs_free` found $\langle\sigma_{FS}\sigma_{PS}^2\rangle$ nonzero/smooth;
  interacting is unexplored.

## What it buys
1. A clean $m_{PS}$ (single-meson) determination from a $2\to1$ amplitude -- independent of the single-$\sigma$
   two-point.
2. The cubic coupling $\langle0|\sigma_{PS}^2|{\rm meson}\rangle$ (how strongly the two-meson interpolator leaks
   to one meson) -- the very contamination that motivated the time-split operators.
3. The one genuine FS-vs-PS test (odd-loop): $T_{FS}\ne T_{PS}$ despite the four-point PS$=$FS.

## Open questions for NM
- Is the target the AMPLITUDE (coupling) itself, or just the $m_{PS}$ cross-check? (sets whether we need
  absolute normalization / $Z$-factors, or just the effmass).
- Include $O_A$ ($=\bar\psi\Phi\tilde\tau\psi$) as a third bilinear for a $\{\sigma_{PS},\sigma_{FS},O_A\}$
  vs $\sigma^2$ comparison, or keep to PS/FS?

## The FS difference and parity (NM 2026-09-09) -- CORRECTED framing
The one place FS $\ne$ PS is this odd (triangle) loop. NM's point (and I had it BACKWARDS at first):

- **Our theory PRESERVES parity.** Parity is an exact symmetry of the target (Nf2, and independent of $N_f$).
  My earlier "survival requires the continuum to break parity / odd-$N_f$ parity anomaly" was WRONG -- there is no
  parity anomaly to invoke here; parity is good.
- **$(1-D_{ov})$ is the explicit parity-BREAKING regulator term.** $\sigma_{PS}=\eta^\dagger\xi+\xi^\dagger\eta$ is
  parity-EVEN; $\sigma_{FS}=\eta^\dagger\xi-\xi^\dagger(1-D_{ov}^\dagger)\eta$ carries the GW factor $(1-D_{ov})$,
  which is exactly the lattice piece responsible for (regulator) parity breaking (same term whose trace gives the
  index). So $\langle\sigma_{PS}^2(t)\,\sigma_{FS}(s)\rangle\ne0$ at finite $a$ is an ARTIFACT of this explicit
  parity-breaking term -- NOT a physical parity-odd correlator (which exact parity forbids).
- **It is testable in the FREE theory.** The nonzeroness is pure kinematics of the $(1-D_{ov})$ term -- no
  dynamics needed. We have `data_free/distill_Nv24` (L1 complete) and `distill_Nv84` (L2); compute
  $T_{FS}-T_{PS}$ there directly and cheaply FIRST.
- **The real open question.** Since parity is exact, $\langle\sigma_{PS}^2\sigma_{FS}\rangle$ is forbidden in the
  continuum and must vanish -- generically $O(a)$, a plain artifact that dies. The interesting question (NM: "not
  sure how much nontrivialness remains in the interacting theory") is whether the INTERACTING theory leaves
  anything beyond the free $O(a)$ artifact -- an anomaly-like protected remnant despite exact parity, or just a
  renormalized $O(a)$ coefficient. Program: (1) FREE-theory value of $T_{FS}-T_{PS}$ (baseline artifact), (2)
  interacting value, (3) continuum scaling at several $L$ -- does it die as $O(a)$ or hold to something finite?

Refs: distillation Peardon 0905.2160; CP mixing 1603.05582; GW index HLN hep-lat/9801021, Luscher hep-lat/9802011.
