# Kernel matrix of the $(2,2)$ state $\{1,1,1,1\}$

**Date:** 2026-09-17. Session "Better understand the {1,1,1,1}". Answers NM's founding question:
what bilinear mesonic operator (kernel) corresponds, via operator-state correspondence, to the
$\{1,1,1,1\}=(2,2)$ level? Companion to `sigma2_single_meson_exclusion_claude.md`,
`sigma_ps_fs_overlap_table_claude.md`, and the free-limit driver `state1111_freelimit_claude.py`.

## The state

$\{1,1,1,1\} = (2,2)$ is the **single-meson** state with **both fermion legs in the $\ell=\tfrac32$
single-particle spinor mode**, contracted to a rotational scalar $J=0$, dimension $\Delta=4$ (continuum).
In the free limit it sits at $a_t m = 2E_{3/2} = 0.556$ (L1) / $\approx 0.69$ (L2), just below the
two-meson threshold $2m_{PS}=4E_{1/2}$ (also $\Delta=4$; the two are degenerate in the continuum, split
only by an $O(a)$ lattice artifact -- gap $0.200$ L1 / $0.096$ L2).

## Single-particle modes: half-integer $\ell$ spinor harmonics

The free propagator $\langle\psi\bar\psi\rangle = D_\text{ov}^{-1}$ on $S^2$ decomposes into **spinor
harmonic modes** labelled by **half-integer** total angular momentum $(\ell,m)$:

$$
\ell = \tfrac12,\ \tfrac32,\ \tfrac52,\ \dots
$$

The half-integer values are forced by the spinor nature of $\psi$: a spin-$\tfrac12$ field is
**anti-periodic** under $\phi\to\phi+2\pi$ (a $2\pi$ rotation acts as $-1$), so its azimuthal modes
$e^{im\phi}$ have half-integer $m$, hence half-integer total $\ell$. This is the direct, complete label --
there is **no** need to split into an integer orbital $L$ plus spin, and **no radial quantum number**
(there is no radial direction on the sphere). The mode index $\ell$ is purely angular.

- $\ell=\tfrac12$ (2 states): the ground meson $m_{PS}=2E_{1/2}$ legs.
- $\ell=\tfrac32$ (4 states): the $(2,2)$ legs.

Energies $E_\ell$ with $m_{PS}=2E_{1/2}$ and $(2,2)=2E_{3/2}$. Going $\ell:\tfrac12\to\tfrac32$ is a higher
**angular** spinor mode, not a radial excitation.

## The kernel

By operator-state correspondence the $(2,2)$ level is created by a local (cylinder-origin) fermion
bilinear

$$
O_{(2,2)}(\tau) = \int_{S^2}\!\!\int_{S^2} \bar\psi(x)\,K(x,y)\,\psi(y) ,
\qquad
K = \Pi_{\ell=3/2} = \sum_{m=-3/2}^{3/2} |{\tfrac32},m\rangle\langle{\tfrac32},m| ,
$$

i.e. the **projector onto the $\ell=\tfrac32$ spinor mode**. This is precisely the $\ell=\tfrac32$
**spectral component of the nonlocal propagator kernel** $D_\text{ov}^{-1}$: in its eigenbasis

$$
D_\text{ov}^{-1} = \sum_\ell \frac{1}{d_\ell}\,\Pi_\ell
\qquad\Longrightarrow\qquad
K_{(2,2)} = \Pi_{\ell=3/2}\ \ (\text{the residue at the }\ell=\tfrac32\text{ mode}).
$$

So the naive guess $\bar\psi\,D_\text{ov}^{-1}\,\psi$ is the *full* propagator kernel summed over all modes;
the $(2,2)$ operator keeps only the $\ell=\tfrac32$ piece on both legs. $K$ is **nonlocal** on $S^2$
because $\Pi_{\ell=3/2}(x,y)=\sum_m u_{3/2,m}(x)\,u_{3/2,m}^\dagger(y)$ connects distinct points -- exactly
the "nonlocal kernel" anticipated at the outset, and the reason the state shows up through the
same-timeslice fermion contraction (diagram A).

## Why the $\ell=0$ (scalar) component exists

The $(2,2)$ meson combines two $\ell=\tfrac32$ spinor modes -- $\psi$ and $\bar\psi$, each in $\ell=\tfrac32$
-- and the scalar is the $J=0$ piece of the product, by the angular-momentum multiplication formula:

$$
\tfrac32 \otimes \tfrac32 = 0 \oplus 1 \oplus 2 \oplus 3 .
$$

The $J=0$ singlet is present, so an $\ell=0$ scalar exists. It is the **rotationally-invariant contraction**
of the two half-integer modes -- the identity/trace on the $\ell=\tfrac32$ multiplet:

$$
K_{\ell=0} = \Pi_{\ell=3/2} = \sum_m |{\tfrac32},m\rangle\langle{\tfrac32},m|\ \ \text{(sum over ALL $m$)} ,
$$

manifestly $SO(3)$-invariant. A single $m$-component, or any partial sum, transforms nontrivially and
carries $J\neq0$; only the **complete** $m$-sum (the identity on the multiplet) realizes the singlet. That
is why the full projector $Q_{\ell=3/2}=\bar\psi\,\Pi_{\ell=3/2}\,\psi$ *is* the rotation-invariant scalar
operator. The other pieces $J=1,2,3$ of $\tfrac32\otimes\tfrac32$ are the non-scalar partners built from the
same two $\ell=\tfrac32$ modes.

## Summary

- **Kernel of the $(2,2)$ state:** $K_{(2,2)} = \Pi_{\ell=3/2}$, the projector onto the $\ell=\tfrac32$
  half-integer spinor mode = the $\ell=\tfrac32$ spectral residue of the nonlocal $D_\text{ov}^{-1}$ kernel.
  Confirms the founding hypothesis that the kernel is a piece of $\bar\psi\,D_\text{ov}^{-1}\,\psi$.
- **Modes are half-integer $\ell$** ($\tfrac12,\tfrac32,\dots$) because spinors are anti-periodic in $\phi$;
  the label is purely angular (no orbital-$L$ split, no radial number). $(2,2)$ is the $\ell=\tfrac32$ mode,
  not a "radial excitation."
- **$\ell=0$ scalar** = the $J=0$ singlet of $\tfrac32\otimes\tfrac32$, realized by the $SO(3)$-invariant
  identity/trace on the shell (the complete $m$-sum in $\Pi_{\ell=3/2}$).
- **Consistency:** $Q_{\ell=3/2}=\bar\psi\,\Pi_{\ell=3/2}\,\psi$ is the clean rotation-invariant scalar
  interpolator, orthogonal to the $\ell=\tfrac12$ ($m_{PS}$) mode; the same state appears as the ground
  ($0.556$ free L1) in the $\sigma^2$ P$_+$ GEVP (`sigma2_3x3_free_gevp_plot_claude.py`, 6x6 state-0), just
  below the two-meson $2m_{PS}=0.756$.
