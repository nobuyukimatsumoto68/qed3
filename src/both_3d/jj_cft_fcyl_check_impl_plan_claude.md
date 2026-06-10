# Implementation plan — check Eq. (4.26), general $n_1\ne n_2$ cylinder current correlator

Equations: `jj_cft_check_equations_claude.md` (qed3int_v2-11). Builds on the validated continuum
propagator and the off-site frame-invariant agreement vs the lattice
(`cont_prop_lattice_compare_claude.md`). Same contraction as the $G_t$ check
(`jj_cft_Gt_samesite_derivation_claude.md`), now at distinct spatial points.

## Goal

Check Eq. (4.26):
$$
f_{\rm cyl}^{ab}(t;\hat n_1,\hat n_2) = C_j\,\Lambda^a_\alpha(\hat n_1)\Lambda^b_\beta(\hat n_2)\,
 e^{-\Delta t}\Big[\delta^{\alpha\beta} - \frac{2 v^\alpha v^\beta}{D}\Big]D^{-\Delta},
\quad v=\hat n_1-e^{-t}\hat n_2,\ \ D=1-2\hat n_{12}e^{-t}+e^{-2t},\ \ \Delta=2,
$$
by contracting the eigenbasis propagator into the vector current,
$f^{ab}=-\mathrm{tr}[\sigma^a G(x_1,x_2)\sigma^b G(x_2,x_1)]$, $x_1=(\hat n_1,t)$, $x_2=(\hat n_2,0)$.

## Key analytic structure (what we actually test, frame-free)

Since $|v|^2 = D$, $\hat v\equiv v/\sqrt D$ is a unit vector and the bracket is the Householder
reflection $H^{\alpha\beta}=\delta^{\alpha\beta}-2\hat v^\alpha\hat v^\beta$ (orthogonal, eigenvalues
$-1,+1,+1$). Thus $T=C_j e^{-\Delta t}D^{-\Delta}H$ and $f=\Lambda_1 T\Lambda_2^{\mathsf T}$ with
$\Lambda_{1,2}\in SO(3)$. Therefore:

1. **Three equal singular values:**
   $$s_1=s_2=s_3 = |C_j|\,e^{-\Delta t}D^{-\Delta} = |C_j|\,e^{-2t}\,(1-2\hat n_{12}e^{-t}+e^{-2t})^{-2}.$$
   This single statement encodes $\Delta=2$, the $D^{-\Delta}$ scaling, AND the reflection tensor.
   Frame-invariant (independent $SO(3)$ on each end), so it is directly comparable to our
   propagator's $f^{ab}$ (whose spinor frame differs from $\Lambda$ by a per-site rotation).
2. **Determinant sign:** $\det f = \det\Lambda_1\det\Lambda_2\,\det T = -\,(|C_j|e^{-2t}D^{-2})^3$
   (reflection has $\det=-1$; $\det\Lambda=+1$).
3. **Temporal component (frame-invariant, fixes orientation):** with the time axis globally $\sigma_3$,
   $$G_t = f^{33} = C_j\,e^{-2t}D^{-2}\big[\hat n_{12} - 2(\hat n_1\!\cdot\hat v)(\hat n_2\!\cdot\hat v)\big],
   \quad \hat n_1\!\cdot\hat v=\frac{1-e^{-t}\hat n_{12}}{\sqrt D},\ \ \hat n_2\!\cdot\hat v=\frac{\hat n_{12}-e^{-t}}{\sqrt D},$$
   reducing to $-C_j e^{-2t}(1-e^{-t})^{-4}$ = Eq. (4.31) at $\hat n_1=\hat n_2$. All quantities depend
   only on $(t,\hat n_{12})$ -> isotropy.

Normalization: with the code's $f=-\mathrm{tr}[\sigma^aG\sigma^bG]$ ($\mathcal N=1$),
$|C_j|=1/(8\pi^2)$ and $C_j=-1/(8\pi^2)$ (fixed by the $\hat n_1=\hat n_2$ case). We CHECK that the
common singular value $/[e^{-2t}D^{-2}]$ is constant $=1/(8\pi^2)$ across all $(t,\hat n_{12})$.

**$G_s$ (general $n_1\ne n_2$) is NOT done component-wise:** $G_s=\sum_a f^{aa}-f^{33}$ uses the full
trace, which is frame-DEPENDENT for $\hat n_1\ne\hat n_2$ (needs the explicit $\Lambda$-frame map).
We rely on the frame-invariants (1)-(3) instead.

## Files

- NEW `jj_cft_fcyl_check_claude.cc` — copies `G_prop`/`xi_mn`(Boost)/`c2_mn`/`cnm2_mn`/`lambda_mn`,
  Pauli $\sigma$, 2x2 helpers; builds the 3x3 $f^{ab}$; computes its singular-value invariants
  ($I_1=\mathrm{tr}\,f^\dagger f$, $I_2$, $I_3=\det f^\dagger f$) without an eigensolver; compares to
  the closed forms. Pure host (GSL + Boost). NO CUDA.
- NEW `run_jj_cft_fcyl_check_claude.sh` — `g++ -x c++` + Boost `-I`, tee to log.

## Chunks

1. **Plumbing.** Copy the propagator functions, $\sigma^a$, `mat2_mul`/`mat2_trace`; add complex 3x3
   helpers (`f^dagger f`, its $I_1,I_2,I_3$, and $\det f$).
2. **Build $f^{ab}$ and invariants.** `fcyl_full(n1,n2,t,nmax) -> cd F[3][3]`; check $\mathrm{Im}\,F\approx0$;
   form $M=F^\dagger F$; $I_1=\mathrm{tr}M$, $I_2=\tfrac12(I_1^2-\mathrm{tr}M^2)$, $I_3=\det M$.
   "All singular values equal $s_0$" iff $I_2=3(I_1/3)^2$ and $I_3=(I_1/3)^3$; then $s_0^2=I_1/3$.
3. **Sweep + compare.** $\hat n_1=(\pi/2,0)$ fixed; $\hat n_2=(\pi/2,\gamma)$ (both equatorial, $z=0$,
   overflow-safe) for $\gamma$ giving $\hat n_{12}=\cos\gamma\in(-1,1)$; $t\in\{0.8,1.2,2.0,3.0\}$.
   Print $\hat n_{12}$, $t$, $D$, $s_0=\sqrt{I_1/3}$, equal-SV residuals, $s_0/[e^{-2t}D^{-2}]$
   (should be const $=1/8\pi^2$), $\mathrm{sign}(\det f)$, $G_t=\mathrm{Re}\,F^{33}$ vs its closed form.
   Add one off-equator pair with matched $\hat n_{12}$ as an isotropy spot-check.

## Numerics

- `nmax=40`, $t\ge0.8$, equatorial $\hat n$ ($z=0$) — the overflow-safe / IR window (the small-$t$ UV
  is not of interest). Single thread; `g++ -x c++ -I../../qfe_mod/include/boost_1_86_0 -lgsl -lgslcblas -lm`.

## Assumptions (same as the $G_t$ check)

Contraction $f^{ab}=-\mathcal N\,\mathrm{tr}[\sigma^aG\sigma^bG]$; temporal axis $=\sigma_3$ globally;
$t=$ cylinder time, $\Delta=2$; overall $C_j$ (sign/magnitude) is convention, checked only via
ratio-constancy / equal-SV / $G_s/G_t$.
