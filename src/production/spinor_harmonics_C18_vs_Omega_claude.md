# Free Dirac spinor harmonics on $S^2$: C.18 vs textbook $\Omega_{jlm}$, and the local-Lorentz map $U(x)$

**Date:** 2026-09-18. Session "Better understand the {1,1,1,1}". Records the exact relationship between the
spinor mode basis used in this codebase (qed3_v2-6.pdf Eq. C.18, `free_wavefunctions_claude.py`) and the
standard "spinor spherical harmonics" $\Omega_{jlm}$ of the literature. Companion to
`state1111_two_factor_kernel_claude.md` (where the $\ell=\tfrac32$ shell projector $P_{3/2}=\sum_a\Xi_a\Xi_a^\dagger$
defines the $(2,2)$ operator). All numbers are free L1 (12-site icosahedron, $N_v=24$ complete basis).

## Two representations of the same eigenspinors

A spinor on $S^2$ has a 2-component index whose meaning depends on the **frame** (vielbein) chosen at each
point. There are two standard choices, and the two "spinor harmonic" conventions correspond to them.

**(A) Global Cartesian frame -- the textbook $\Omega_{jlm}$.** Refer every spinor to the single fixed frame of
the embedding $\mathbb{R}^3$ ($\hat x,\hat y,\hat z$), the same at every point. Then $\chi_\pm$ (eigenstates of
$\sigma_z$ along the constant global $\hat z$) are position-independent, and the spinor harmonic is the
Clebsch-Gordan coupling of an **ordinary (integer-$\ell$, periodic) $Y_{\ell m}$** to $\chi_\sigma$:
$$
\Omega_{j,\,\ell=j\mp1/2,\,m}(\theta,\phi)=\sum_{\sigma=\pm1/2}\langle \ell,\,m-\sigma;\ \tfrac12,\sigma\,|\,j,m\rangle\;Y_{\ell,\,m-\sigma}(\theta,\phi)\,\chi_\sigma .
$$
The $j=\tfrac32$ shell is spanned by $\{\Omega_{3/2,1,m},\Omega_{3/2,2,m}\}_{m=-3/2}^{3/2}$ (both $\ell=1,2$, 8
functions). Standard refs: Edmonds (1957); Varshalovich–Moskalev–Khersonskii (1988).

**(B) Local tangent frame -- C.18 (this codebase).** Refer the spinor to the *moving* orthonormal frame
$(e_\theta,e_\phi,e_r)$ intrinsic to $S^2$. The curved-space Dirac operator lives here (it needs the vielbein +
spin connection). Its eigenspinors are the spin-weighted / Jacobi form, `free_wavefunctions_claude.psi`:
$$
\psi_{m,n,\iota_3}\propto e^{im\phi}\begin{pmatrix}\xi_{|m|,n}(\iota_m z)\\[2pt]\iota_3\,\iota_m\, i\,(-1)^n\,\xi_{|m|,n}(-\iota_m z)\end{pmatrix},\quad
\xi_{|m|,n}(z)=(1-z)^{\frac{|m|-1/2}{2}}(1+z)^{-\frac{|m|+1/2}{2}}P^{(|m|-1/2,-|m|-1/2)}_{n+|m|+1/2}(z),
$$
with $z=\cos\theta$, $\iota_m=\mathrm{sign}(m)$, half-integer $m$. The frame rotates by $2\pi$ around the
azimuth, so the spinor is **anti-periodic** in $\phi$ (the half-integer $e^{im\phi}$) -- this is the
spin-weighted spherical harmonic / Wu–Yang $q=\pm\tfrac12$ monopole-harmonic tradition. Standard ref for the
Dirac eigenspinors on $S^n$: Camporesi & Higuchi, J. Geom. Phys. **20** (1996) 1, arXiv:gr-qc/9505009.

## Finding: same eigenspace, different frame (local Lorentz rotation)

The two shell projectors $P^{\rm C.18}_{3/2}=\sum_a\Xi_a\Xi_a^\dagger$ and
$P^{\Omega}_{3/2}=\sum\Omega\,\Omega^\dagger$ (each the rank-8 orthogonal projector onto its $j=\tfrac32$ modes)
are **different projectors** ($\|P^{\rm C.18}-P^\Omega\|/\|P\|\approx1.0$) but describe the **same physical
eigenspace in two spin frames**, related by a per-site local Lorentz rotation:
$$
\boxed{\;P^{\rm C.18}(x,z)=U(x)\,P^{\Omega}(x,z)\,U(z)^\dagger,\qquad \Xi(x)=U(x)\,\Omega(x)\ (\text{mod a within-shell rotation}),\quad U(x)\in SU(2).\;}
$$

**Numerical evidence (all free L1):**
- **Frame-invariant density** $\mathrm{Tr}_\text{spin}P(x,x)=0.6667$ constant for both, matching to $10^{-16}$
  ($P(x,x)\propto I_2$, spin-isotropic, as a full multiplet must be).
- **2-point invariant** $\mathrm{Tr}[P(x,z)P(z,x)]$ and **3-point invariant** $\mathrm{Tr}[P(x,y)P(y,z)P(z,x)]$
  match between C.18 and $\Omega$ to $\sim10^{-16}$ (frame-invariant, so equal iff same eigenspace up to a
  per-site rotation).
- **$U(x)$ constructed explicitly** by Procrustes synchronization (init $U=I$; iterate
  $U(x)=\mathrm{polar}\!\sum_z P^{\rm C.18}(x,z)U(z)P^{\Omega}(z,x)$): converged in 57 iters with
  $$
  \max_{x,z}\big\|\,U(x)\,P^{\Omega}(x,z)\,U(z)^\dagger-P^{\rm C.18}(x,z)\,\big\|=3.4\times10^{-14},
  $$
  each $U(x)$ unitary with $\det=1$ ($SU(2)$) to $10^{-16}$. (A naive single-reference-column extraction fails
  because the $2\times2$ blocks $P(x,x_0)$ are rank-deficient; the synchronization fixes that.)

So it is a genuine local Lorentz (spin-frame) relation of the **spinors**, not merely a projector or density
coincidence.

## What $U(x)$ is

$U(x)$ is the **vielbein / spin-connection rotation** ($D^{1/2}$) from the global Cartesian frame to the local
$S^2$ tangent frame at $(\theta,\phi)$,
$$
U(\theta,\phi)=e^{-i\phi\,\sigma_3/2}\,e^{-i\theta\,\sigma_2/2}\,e^{+i\phi\,\sigma_3/2}\quad(\text{up to a global }SU(2)\text{ gauge}),
$$
the $e^{\pm i\phi\sigma_3/2}$ being the azimuthal frame twist that converts the periodic integer-$\ell$ $Y_{\ell m}$
of $\Omega$ into the anti-periodic half-integer-$m$ spin-weighted form of C.18. (The exact Euler convention is
fixable from the axis–angle of the numerical $U(x)$; the constructed $U(x)$ carries a residual global-$SU(2)$
gauge.)

## Consequences / usage

- **They are the same physics.** Any shell projector, density, or meson correlator built from C.18 equals the
  one built from $\Omega$ (and from the stored lattice eigenbasis $V$ -- see below): all are the frame-invariant
  $j=\tfrac32$ shell. The $(2,2)$ operator $\bar\psi\,P_{3/2}\,\psi$ is frame-independent.
- **Lattice caveat.** Both C.18 and $\Omega$ are *continuum* modes; on the coarse 12-site lattice the free
  $D_\text{ov}$ has only icosahedral symmetry, so a single continuum-shell projector is ground-contaminated in
  the meson correlator (state1 obstructed at L1, clean at L2). The lattice's OWN $\ell=\tfrac32$ shell -- the
  8-fold eigenspace of $M=\tau(t_0,t_0)-\tfrac12$ (stored perambulator) or the $\lambda=2$ columns of the stored
  $/V$ by $/\text{evals}$ -- gives the clean $2E_{3/2}=0.556$ without any $\Xi$/$\Omega$ object or extra solve
  (see `state1111_two_factor_kernel_claude.md`).

## References
- R. Camporesi, A. Higuchi, *"On the eigenfunctions of the Dirac operator on spheres and real hyperbolic
  spaces,"* J. Geom. Phys. **20** (1996) 1–18, arXiv:gr-qc/9505009. (Dirac eigenspinors on $S^n$ -- the C.18 /
  spin-weighted form.)
- A. R. Edmonds, *Angular Momentum in Quantum Mechanics*, Princeton (1957); D. A. Varshalovich, A. N. Moskalev,
  V. K. Khersonskii, *Quantum Theory of Angular Momentum*, World Scientific (1988). (CG-coupled $\Omega_{jlm}$.)
- T. T. Wu, C. N. Yang, Nucl. Phys. **B107** (1976) 365. (Monopole harmonics; reduce to spinor harmonics at
  $q\to0$; the tradition used in e.g. arXiv:2112.02106.)

Driver / checks: `state1111_operator_verify_claude.py` (C.18 shell modes, `build_shell_Xi`, `rot_and_su2`),
inline scripts for the $\Omega$ construction, the frame invariants, and the Procrustes $U(x)$.
