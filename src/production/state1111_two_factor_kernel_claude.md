# The $\{1,1,1,1\}=(2,2)$ operator: two-factor $\ell=\tfrac32$ kernel, verified free

**Date:** 2026-09-18. Session "Better understand the {1,1,1,1}". Confirms NM's precise operator picture for
the $(2,2)$ single-meson state in the free limit. Driver: `state1111_operator_verify_claude.py`. Supersedes
the kernel description in `state1111_kernel_matrix_claude.md` (which had the correct shell but the wrong
"single projector" form).

## The operator (NM)

The $\{1,1,1,1\}=(2,2)$ one-meson operator is an **equal-time nonlocal bilinear** whose kernel is the
$\ell=\tfrac32$ propagator **squared** — one factor dressing each fermion field:
$$
O_{\{2,2\}}(t) = P_{\ell=0}\Big[\ \bar\psi(t,x)\,\Xi_{3/2}(x,y)\,\Xi_{3/2}(y,z)\,\psi(t,z)\ \Big],
\qquad \Xi_{\ell} := D^{-1}\big|_{\ell},
$$
with the spatial indices $x,y,z$ summed on the timeslice $t$. Here $\Xi_\ell$ is the **2D spinor propagator
restricted to the $\ell$-shell** — the C.18 spinor basis of the free-limit paper (qed3_v2-6.pdf Eq C.18),
built in code by `free_wavefunctions_claude.psi(m,n,i3)`. Key points:

- **$\Xi_{3/2}$ is not a projector.** It is the propagator on the shell and carries the Dirac/energy-sign
  structure ($\pm\lambda$ for the $\iota_3=\pm$ particle/antiparticle spatial modes),
  $\Xi_{3/2}=\tfrac{1}{\lambda}(P_+ - P_-)$.
- **Two factors, one per fermion field**, put **both** internal lines in $\ell=\tfrac32$ ("internally
  rotating" — the on-timeslice $x\to y\to z$ propagation carries angular momentum $\tfrac32$). On a single
  degenerate shell $\Xi_{3/2}^2=\tfrac{1}{\lambda^2}(P_++P_-)\propto$ the shell density
  $P_{3/2}=\sum_{m,\iota_3}\Xi_a\Xi_a^\dagger$, so the two-factor kernel reduces to the shell projector.
- **The physical $0^{++}$ is the scalar ($\ell=0$) projection** of the "complex" (tensor) operator: the full
  $\bar\psi\,\Xi_{3/2}^2\,\psi$ carries the whole $\tfrac32\otimes\tfrac32=0\oplus1\oplus2\oplus3$ tower, and
  the $(2,2)$ is its $J=0$ singlet (the higher $J$ are the non-scalar partners of the same shell).

Energy: both constituents in shell $\ell$ $\Rightarrow$ meson energy $2E_\ell$. So $\ell=\tfrac12$ gives
$m_{PS}=2E_{1/2}$ and $\ell=\tfrac32$ gives $(2,2)=2E_{3/2}$.

## The outer $\ell=0$ projector (Clebsch-Gordan), $m$-resolved

Promoting the shell magnetic label, $\Xi^{m}_{3/2}$ = the individual $j=\tfrac32$ multiplet member. The two
dressed legs each carry $j=\tfrac32$, so the $\ell=0$ scalar is the $J=0$ Clebsch-Gordan singlet of
$\tfrac32\otimes\tfrac32=0\oplus1\oplus2\oplus3$ — **couple $m$ with $-m$**, weighted by the CG:
$$
|0,0\rangle=\sum_{m=-3/2}^{3/2}\langle\tfrac32\,m;\tfrac32\,{-m}|00\rangle\,|\tfrac32,m\rangle\otimes|\tfrac32,-m\rangle,
\qquad \langle\tfrac32\,m;\tfrac32\,{-m}|00\rangle=\frac{(-1)^{3/2-m}}{2},
$$
$$
|0,0\rangle=\tfrac12\big(|\tfrac32,\tfrac32\rangle|\tfrac32,-\tfrac32\rangle-|\tfrac32,\tfrac12\rangle|\tfrac32,-\tfrac12\rangle
+|\tfrac32,-\tfrac12\rangle|\tfrac32,\tfrac12\rangle-|\tfrac32,-\tfrac32\rangle|\tfrac32,\tfrac32\rangle\big).
$$
So the kernel is NM's $\sum_m\Xi^m_{3/2}\,\Xi^{-m}_{3/2}$, completed with the CG sign $(-1)^{3/2-m}$ and the
charge-conjugation that turns the $\bar\psi$ (conjugate rep) leg into a $|{-m}\rangle$ ket. The C.18 spinors
satisfy the exact time-reversal relation (verified to $10^{-16}$ on the grid):
$$
\Xi^{-m}_{3/2}(x)=i\,(-1)^{3/2-m}\,(i\sigma_2)\,\big[\Xi^{m}_{3/2}(x)\big]^{*},
$$
and with it the $m\leftrightarrow-m$ CG kernel **equals the shell projector** — verified as a matrix identity
on the full 8-mode $\lambda=2$ shell (residual $1.5\times10^{-16}$):
$$
P_{3/2}\;=\;\sum_{\iota_3=\pm}(i\,\iota_3)\sum_{m}(-1)^{3/2-m}\,\Xi^{m,\iota_3}_{3/2}\otimes\big(i\sigma_2\,\Xi^{-m,\iota_3}_{3/2}\big).
$$
Substituting the conjugation to eliminate $\Xi^{-m}$ (using $i\sigma_2\,\Xi^{-m,\iota_3}=-(i\iota_3)(-1)^{3/2-m}(\Xi^{m,\iota_3})^{*}$)
the CG sign $(-1)^{3/2-m}$, the charge-conjugation $i\sigma_2$, and the $i\iota_3$ phase **all cancel**, leaving simply
$$
P_{3/2}\;=\;\sum_{m,\iota_3}\,\Xi^{m,\iota_3}_{3/2}\,\big(\Xi^{m,\iota_3}_{3/2}\big)^{*}
\;=\;\sum_a \Xi_a\,\Xi_a^{\dagger},
$$
the plain $m$-diagonal shell projector — the form used in the driver. So writing the conjugate leg as $\Xi^{*}$
absorbs the whole CG apparatus: the $m\leftrightarrow-m$ coupling with its $(-1)^{3/2-m}$ sign and $i\sigma_2$ is
*identically* the plain $\sum_a\Xi_a\Xi_a^{\dagger}$. This is why $\bar\psi\,P_{3/2}\,\psi$ was already the $\ell=0$ scalar: the conjugate-rep trace $\sum_m\bar\psi_m\psi_m$
**is** the $m\leftrightarrow-m$ CG contraction, the $(-1)^{3/2-m}$ and the $i\sigma_2$ being exactly the phases
that relate $\langle\Xi^m|$ to $|\Xi^{-m}\rangle$. (The $J=1,2,3$ pieces of $\tfrac32\otimes\tfrac32$ are the
non-scalar partners, obtained with the corresponding higher-$J$ CG weights instead.)

## Construction on the lattice

A meson two-point with a spatial kernel $K$ reduces to the standard distillation trace,
$$
\langle O(t)\,O(0)\rangle = -\mathrm{Tr}\!\big[\Phi_K(t)\,\tau(t,s)\,\Phi_K(s)\,\tau(s,t)\big],
\qquad \Phi_K = V^\dagger K V ,
$$
so for $K=P_\ell$ the vertex is $\Phi_\ell(t)=O O^\dagger$ with $O[c,a]=\langle V_c(t)\,|\,\Xi_a\rangle$
(area-weighted overlap of the distillation vectors with the C.18 shell modes). The L1 icosahedron has two
vertices exactly at the poles, where the C.18 spinor is coordinate-singular; since $P_\ell$ is
rotation-invariant, the modes are evaluated in a **rotated frame** (with the SU(2) spinor rotation $U_R$) so
no site hits a pole — verified by mode-overlap orthonormality to $10^{-16}$.

## What we found

1. **Continuum shells are clean.** The C.18 $\ell=\tfrac12,\tfrac32,\tfrac52$ modes are exactly orthonormal
   on the lattice grid ($\|$off-diag$\|\sim10^{-16}$; $\lambda=1,2,3$ together form a full orthonormal basis).
2. **But the free *lattice* $D_\text{ov}$ mixes them.** The lattice keeps only **icosahedral** symmetry, not
   full $SO(3)$, so $\tau$ is not shell-diagonal in the continuum basis (off-shell block $\sim$ the diagonal
   at L1). A *single* $P_\ell$ operator is therefore ground-contaminated: both $P_{1/2}$ and $P_{3/2}$
   two-points decay to the ground $m_{PS}$ at large $t$. This is a coarse-lattice artifact, not a failure of
   the picture (the clean-shell statement is continuum).
3. **A GEVP of shell operators separates the states.** Basis $\{\,\bar\psi P_\ell\psi : \ell=\tfrac12,\tfrac32,(\tfrac52)\,\}$,
   GEVP at $T_0=2$:

   | | state0 (predict $m_{PS}$) | state1 (predict $(2,2)=2E_{3/2}$) |
   |---|---|---|
   | **L1** (12 sites) | $0.378$ (plateau $t8$-$14$) $=m_{PS}$ | $\approx0.43$ — **obstructed** (transits $\sim0.5$ at $t4$-$5$ then shell-mix contamination) |
   | **L2** (42 sites) | $0.40\approx m_{PS}=0.393$ | $\mathbf{0.693}$ (plateau $t6$-$9$) $=(2,2)=0.690$ |

   **state0 is $m_{PS}$ at both $L$** (both legs $\ell=\tfrac12$). **state1 lands on $2E_{3/2}=0.69$ at L2**
   (both legs $\ell=\tfrac32$); at coarse L1 the icosahedral shell-mixing obstructs it (it only transits the
   $(2,2)$ region). The isolation **sharpens with refinement** — exactly as expected for a continuum
   shell picture realized on a discrete lattice.

## Conclusion

The picture is confirmed: **$\{1,1,1,1\}=(2,2)$ is a single meson with both constituents in the
$\ell=\tfrac32$ shell, energy $2E_{3/2}$** (L2: $0.693$ vs predicted $0.690$). Its operator is the equal-time
nonlocal bilinear $\bar\psi\,\Xi_{3/2}^2\,\psi$ with $\Xi_{3/2}=D^{-1}|_{\ell=3/2}$ the C.18 spinor
propagator on the shell (two factors, one per leg, "internally rotating"), and the physical $0^{++}$ is its
$\ell=0$ scalar projection — the $J=0$ singlet of $\tfrac32\otimes\tfrac32$. On the lattice the clean shell
is only approximate (icosahedral $\subset SO(3)$), so the state is isolated variationally (shell GEVP) and
sharpens as $L$ grows; the free L2 result confirms it quantitatively.

Figures: `figs/state1111_operator_verify_L1_claude.png`, `figs/state1111_operator_verify_L2_claude.png`.
