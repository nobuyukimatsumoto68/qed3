# Continuum all-to-all fermion propagator from the eigenbasis (Eq. C.28) — implementation plan

## Physics / goal

Evaluate the **exact continuum** free fermion propagator on $S^2\times\mathbb{R}$ in the
eigenbasis, Eq. (C.28) of `qed3_v2-6.pdf` (label `eq:Gf` in `qed3_v2(2)/main.tex`):

$$
G(x,x')=\sum_K \frac{\psi_K(x)\,\psi_K^\dagger(x')}{i\,\Lambda_K\,c_K^2},
\qquad K=(k,m,n,\iota_3),
$$

at the **lattice site coordinates** of the simplicial $S^2$ for $L=1,2,4$, to serve as the
analytic reference for the other agent's exact lattice all-to-all (LU) propagator.

Reference for the eigenvalue construction: G. 't Hooft / Abrikosov-style Dirac-on-$S^2$,
`Abrikosov:2002jr` (cited in the appendix). The temporal specialization Eq. (C.29) is the
validation anchor.

### Spectrum and eigenfunctions (from App. C)

$$
\Lambda_{|k|,|m|,n}=\sqrt{k^2+M^2},\qquad M\equiv n+|m|+\tfrac12,
\qquad m\in\mathbb{Z}+\tfrac12,\ n\in\mathbb{Z}_{\ge0},\ \iota_3=\pm1,
$$
$$
\psi_{k,m,n,\iota_3}(\theta,\phi,t)=\frac{e^{\iota_3 i k t}}{\sqrt{2\pi}}\,
\frac{e^{i m\phi}}{\sqrt{2\pi}}\,\chi^{(\mathrm{sgn}\,k)}_{k,m,n,\iota_3}(\theta),
$$
$$
\chi^{(k>0)}=\begin{pmatrix}\cos\frac{\gamma}{2}\,\xi_{|m|,n}(\iota_m z)\\[2pt]
\iota_3\iota_m\,i(-1)^n\sin\frac{\gamma}{2}\,\xi_{|m|,n}(-\iota_m z)\end{pmatrix},
\quad
\chi^{(k<0)}=\begin{pmatrix}\sin\frac{\gamma}{2}\,\xi_{|m|,n}(\iota_m z)\\[2pt]
\iota_3\iota_m\,i(-1)^n\cos\frac{\gamma}{2}\,\xi_{|m|,n}(-\iota_m z)\end{pmatrix},
$$
with $z=\cos\theta$, $\iota_m=\mathrm{sign}(m)$, $\cos\gamma=|k|/\Lambda$, so
$\cos^2\frac\gamma2=\tfrac12(1+|k|/\Lambda)$, $\sin^2\frac\gamma2=\tfrac12(1-|k|/\Lambda)$,
$\cos\frac\gamma2\sin\frac\gamma2=\tfrac12\,M/\Lambda$.

Angular function (Jacobi polynomial form, App. C):
$$
\xi_{|m|,n}(z)=(1-z)^{\frac12(|m|-\frac12)}(1+z)^{-\frac12(|m|+\frac12)}
P^{(|m|-\frac12,\,-|m|-\frac12)}_{\,n+|m|+\frac12}(z),
$$
$$
c_{|m|,n}^2=\frac{1}{2n+2|m|+1}\,
\frac{n!\,(n+2|m|)!}{(n+|m|-\frac12)!\,(n+|m|+\frac12)!}.
$$
Edge values (Eq. C.15 `eq:xi_deltam`): $\xi_{|m|,n}(1)=\delta_{|m|,1/2}/\sqrt2$, $\xi_{|m|,n}(-1)=0$.

### Reduction over the continuous label $k$

$\sum_K\to\sum_{\iota_3=\pm1}\sum_{m}\sum_{n}\int_{-\infty}^{\infty}\frac{dk}{2\pi}$.
Every $k$-dependent factor in $\psi_K\psi_K^\dagger/(i\Lambda)$ is $e^{ik\,\tau}$ ($\tau=t-t'$)
times one of $\{1/\Lambda,\ |k|/\Lambda^2,\ M/\Lambda^2\}$. The $\iota_3$ sum symmetrizes
$\tau\to\pm\tau$. The three building-block $k$-integrals are standard:
$$
\int_{-\infty}^{\infty}\!\frac{dk}{2\pi}\,\frac{e^{ik\tau}}{\Lambda}
=\frac{1}{\pi}K_0(M|\tau|),\qquad
\int_{-\infty}^{\infty}\!\frac{dk}{2\pi}\,\frac{M\,e^{ik\tau}}{\Lambda^2}
=\tfrac12\,e^{-M|\tau|},
$$
$$
\int_{-\infty}^{\infty}\!\frac{dk}{2\pi}\,\frac{|k|\,e^{ik\tau}}{\Lambda^2}
=-\frac{1}{\pi}\big[\,\cos(M\tau)\,\mathrm{Ci}(M|\tau|)+\sin(M|\tau|)\,\mathrm{si}(M|\tau|)\,\big]
\ \ (\text{the }\mathrm{Ci}\text{-type term}).
$$
All available in GSL: `gsl_sf_bessel_K0`, `exp`, `gsl_sf_Ci`, `gsl_sf_Si`.

**Implementation choice for the $k$-integral:** do it by **GSL adaptive quadrature**
(`gsl_integration_qag`/`qagi`) of the single smooth 1D integrand, uniformly for all terms.
Rationale: avoids branch/special-function bookkeeping of the $|k|/\Lambda^2$ piece and the
$\iota_3$/$k$-sign casework; the integrand decays and is smooth for $\tau\ne0$. The
closed-form Bessel/Ci expressions above are kept in comments as a cross-check and as a later
optimization. (Validate both against Eq. C.29.)

After the $k$-integral, $G$ is a **double sum over $n\ge0$ and half-integer $m$** of
[angular factors $\xi$ at $\theta,\theta'$] $\times$ [$e^{im(\phi-\phi')}$] $\times$
[$1/c^2$] $\times$ [the $\tau$-kernel above]. Convergence is exponential in $M$ for
$\tau\ne0$ (factor $e^{-M|\tau|}$ / $K_0(M|\tau|)$); truncate at $n_{\max},|m|_{\max}$ with a
relative-tail stopping criterion.

### Validation anchor (temporal, Eq. C.29)

For $x=(0,0,t)$, $x'=(0,0,0)$ only $|m|=1/2$ survives ($\xi(1)=\delta_{|m|,1/2}/\sqrt2$) and
$$
G(t)=\mathrm{sign}(t)\,\sigma_3\sum_{n\ge0}\frac{n+1}{4\pi}\,e^{-(n+1)|t|}.
$$
The general evaluator at $\theta=\theta'=0$ must reproduce this (and the $\Delta=1$ tower).

## Scope decisions (from user)

- **Separation-reduced:** exploit $t$-translation (depends on $\tau=t-t'$ only) and the fact
  that $G$ depends on $(\theta,\theta',\phi-\phi',\tau)$. Compute over inequivalent
  unordered site pairs $\times$ $\tau\in\{0,\dots,N_t-1\}$, i.e. $\binom{N_s+1}{2}\,N_t$
  evaluations instead of $(N_s N_t)^2$; expand to the full all-to-all matrix on output.
- **Geometry from disk:** read site Cartesian coords from
  `geometry/data/pts_n{1,2,4}.dat` (unit sphere; 12/42/162 lines), convert to
  $(\theta,\phi)$ via $\theta=\arccos z$, $\phi=\mathrm{atan2}(y,x)$.

## Resolved decisions

1. **$t$ normalization / range — RESOLVED.** Use the `hmc_claude.cu` values:
   $a_t = 0.2$ (`hmc_claude.cu:201`), $N_t = 128$ (`hmc_claude.cu:69`). Continuum
   $t = a_t\times(\text{time index})$, index $0\dots N_t-1$, sphere radius $=1$ (sites read
   on the unit sphere). $T = a_t N_t = 25.6$. CLI-overridable; $L$ (refinement) is a
   parameter, $a_t,N_t$ default to the above.
2. **Output — RESOLVED: match the other project's `Dinv.h5` schema exactly** (user). From
   `jj_propagator_exact_claude.cu` (write block ~L362-378), one file
   `data_<esnid>/prop_exact_L<L>/Dinv.<k>.h5` per $L$ with datasets:
   - `Dm_inv/real`, `Dm_inv/imag` — **row-major length-$N^2$**, index $=\text{row}\,N+\text{col}$,
     value $=(D_m^{-1})_{\text{row,col}}$. We fill this with the continuum $G$ (Sec. index map).
   - `N`, `N_REFINE`, `Nt`, `n_sites` (ints), `parity` (int), `mP/{real,imag}`, `k` (int),
     `complete` (int sentinel, written LAST), and `Dtil_inv/{real,imag}` if `parity`.
   - atomic write: `.h5.tmp` then `std::filesystem::rename`.
   Build flags (from `run_reweighting_R_claude.sh`):
   `-I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/`,
   link `-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -lhdf5 -lgsl -lgslcblas -lm`.
   So the **same** `jj_contract_exact_claude.cu` / notebooks read our continuum file identically.
3. **Code location — `src/both_3d/`**, standalone `cont_prop_eigbasis_claude.cpp`
   (CPU + GSL + HighFive), built/run via `run_cont_prop_claude.sh`.
4. **Global DOF index map — RESOLVED** (`includes/dirac_ext.h:65`,
   `is = Nx*s + NS*ix + a`): linear index
   $$ r = N_x\,s + N_S\,i_x + a,\qquad N_x=N_S N_{\rm sites},\ N_S=2, $$
   i.e. **time $s$ slowest, site $i_x$ middle, spinor $a\in\{0,1\}$ fastest**. Decode:
   $a=r\bmod2$, $i_x=(r/2)\bmod N_{\rm sites}$, $s=r/(2N_{\rm sites})$. Site index $i_x$ uses
   the SAME ordering as `geometry/data/pts_n<L>.dat` lines (the base lattice site order).
   We fill `Dm_inv[r N + c]` $=G^{ac}(x_{i_x,s},x_{i_y,s'})$ for the decoded pairs.
5. **nmax — RESOLVED: $n_{\max}=40$** (user); $|m|$ truncation likewise bounded (the $\xi$
   factor vanishes unless $|m|\le n+\tfrac12$ range contributes; tail set by $e^{-M|\tau|}$).
6. **Time separations — RESOLVED: only $t\ne t'$** (user): we ignore the contact term and
   only evaluate correlators between two DIFFERENT time slices. Then $\tau=t-t'\ne0$ always,
   so the $(n,m)$ sum carries full $e^{-M|\tau|}$ suppression — uniformly convergent, no
   coincident/short-distance subtlety. (Same-site, same-time diagonal blocks: leave unset/zero.)
7. **Overall normalization — RESOLVED: irrelevant** (user). We drop the $\bar a_s a_t$ and
   $1/4\pi$ prefactors; only the functional form / relative structure matters. The
   normalization subsection below is kept for reference only.

## Integration with the other agent (exact lattice all-to-all)

The other agent's plan is `jj_exact_freefield_impl_plan_claude.md`. What it does and how we dock:

- **Program** `jj_propagator_exact_claude.cu` (sibling of `reweighting_R_claude.cu`): dense build
  of $D_{\rm ov}$ by `Op.from_cpu<N>(Dxi,e_i)` (one column per unit vector), then cuSOLVER
  `Zgetrf`/`Zgetrs` LU inverse -> saves the **lattice all-to-all** $[D_{\rm lat}^{-1}]_{xx'}$ to
  `data_<esnid>/prop_exact_L<L>/Dinv.<k>.h5` (row-major, complex, atomic write). Free field
  ($U=1$) => single config, built once. Index size $N=n_{\rm sites}\times N_t\times 2$
  (L1 3072 / L2 10752 / L4 41472). L1/L2 fit one GPU; L4 (27.5 GB) needs A100-80 + host trims.
- **Program** `jj_contract_exact_claude.cu`: loads `Dinv`, contracts into the **gauge-invariant**
  current correlators (tp/sp/ylm, vector + axial), validates Eqs. (4.28)/(4.31)/(4.34) and the
  parameter-free ratio $G_s=-(D-1)G_t$.

**Our role:** Eq. (C.28) is the *continuum* analytic reference these lattice numbers converge to.
Two clean comparison levels (both gauge invariant in free field):
1. **Propagator level** (the `prop_ov` / `comp_ferm_2d` checks): our $G(x,x')$ vs their
   $G_{\rm lat}=\frac{1}{\bar a_s a_t}[D_{\rm lat}^{-1}]_{xx'}$.
2. **Correlator level**: feed our continuum $G$ through the *same* gauge-invariant contractions
   to predict (4.28)/(4.31)/(tower) directly — a continuum cross-check of their contraction stage.

### Normalization (paper Sec. `sec:dirac_propagators`)
$$
G_{\rm lat}(x,x') = \frac{1}{\bar a_s\,a_t}\,[D_{\rm lat}^{-1}]_{xx'}\ \xrightarrow{\text{cont}}\ G(x,x')\ \text{(Eq. C.28)}.
$$
So to compare against their **raw** `Dinv.h5`, scale: $[D_{\rm lat}^{-1}]\approx \bar a_s a_t\,G_{\rm cont}$.
- $a_t = 0.2$ (`hmc_claude.cu:201`).
- $\bar a_s = $ `base.mean_ell` (mean simplicial edge length; **per $L$**, printed as `# alat`
  at `hmc_claude.cu:239`; constraint `sqrt(3) abar_s/a_t >= 4/sqrt(3)`, `hmc_claude.cu:202`).
  Read it from the geometry/lattice setup at the matching $L$ (do NOT hardcode).

### Reference code: `dirac_inverse/integrate.cc` (user's earlier $S^2$ implementation)
Same eigenbasis, **identical $\xi$ convention** (validated against ours):
- `xi(mpH,n,z)`, `mpH = |m|+1/2` (integer): factor $(1-z)^{(mpH-1)/2}(1+z)^{-mpH/2}$ times
  `boost::math::jacobi(n+mpH, mpH-1, -mpH, z)` = our $P^{(\alpha,\beta)}_j$ with
  $\alpha=mm-\tfrac12,\ \beta=-(mm+\tfrac12),\ j=n+mm+\tfrac12$, $mm=|m|=mpH-\tfrac12$.
- `Cnm(mpH,n) = sqrt( 4*pi*c2_{|m|,n} )` — their normalization carries an extra $4\pi$ (the
  overall factor we drop). Our `c2_mn` equals `Cnm^2/(4*pi)`.
- Commented `psiP/psiM`: the $2$-spinor $(\xi(z),\,i(-1)^n\xi(-z))$ with phase
  $e^{i\,\mathrm{sign}(m)(mpH-1/2)\phi}$ over `Cnm` — the App. C $\psi_{m,n,\iota_3}$.
- `main` line ~138 IS the $S^2$ CFT check: $\sum_n(-1)^n\xi(1,n,-z)/(\sqrt2\,\pi)$ vs
  $1/(4\pi\sin(\theta/2))$ ($(1,2)$ comp). Reuse verbatim for our Chunk 3 $S^2$ test.
- Boost is available under `dirac_inverse/` (uses `boost/math/special_functions/jacobi.hpp`)
  as an INDEPENDENT Jacobi cross-check if wanted.

### Analytic anchors for validation (besides C.29)
- **Temporal** (Eq. `eq:exact_Gf`, the $(1,1)$ component, `prop_ov` figure):
  $G(t)=\sigma_3\,\mathrm{sign}(t)\,\frac1{4\pi}\sum_{n\ge0}(n+1)e^{-(n+1)|t|}$.
- **Pure $S^2$ CFT** (Eqs. `eq:psipsi_S2`/`eq:psipsi_S2_ev`, the $(1,2)$ component,
  `comp_ferm_2d` figure): drop the time direction; the eigenbasis sum
  $\sigma_1\sum_{n\ge0}\frac{(-1)^n}{\pi\sqrt2}\xi_{1/2,n}(-z)$ must approach the closed form
  $\dfrac{\sigma_1}{4\pi\sin(\theta/2)}$ as $n_{\max}\to\infty$ ($n_{\max}=50,200$ in the paper).
  This directly validates our $\xi$ primitive (Chunk 1) against a known function.

## Time boundary condition — RESOLVED (user)

No antiperiodic folding. Use the literal $\tau = a_t\,(s-s')$ and evaluate $G$ for both signs of
$t$ out to $|\tau| = a_t N_t$. Periodic images are dropped (negligible, $e^{-\lambda T/2}$). The
$s=s'$ ($\tau=0$) diagonal-in-time blocks are left zero (contact term, not wanted).

## Chunks 1-4 — DONE and validated (`cont_prop_eigbasis_claude.cpp` + `run_cont_prop_claude.sh`)

- C1 `xi_mn`, `c2_mn`: edge values 6e-16, orthonormality $\int\xi\xi=c^2$ to 3e-6 ($n\le40$).
- C2 `G_prop` full $2\times2$ (k-integral done analytically -> pure exponentials
  $e^{-\lambda_{|m|,n}|\tau|}$, prefactor $1/(4\pi c^2)=1/C_{nm}^2$): reproduces (C.29) to 8e-16,
  exact $\sigma_3$.
- C3 $S^2$ CFT check 4e-6 at $\theta=\pi/2$; 3D self-convergence in $n_{\max}$ to 1e-9
  ($\tau\ge0.5$), nearest-slice $\tau=0.2$ truncates at $\sim e^{-8}$ as predicted.
- C4 `build_all_to_all`: reads `pts_n<L>.dat`, precomputes $A,B$ spinor tables, separation-reduced
  scatter into the row-major $N^2$ matrix (both $\tau$ signs, $\tau=0$ left zero), writes
  `cont_prop_L<L>/Dinv.<k>.h5` in the EXACT `Dinv.h5` schema. Validated: L=1 (0.9 s, 0.15 GB) and
  L=2 (10 s, 1.85 GB) reproduce the north-pole same-site (C.29) block exactly, negative-$t$ flips
  the diagonal (odd), $\tau=0$ blocks zero. Run `bash run_cont_prop_claude.sh [1 2 4]`; CLI
  `--L --nt --at --nmax --k --geo --out`. L=4 (N=41472, 27.5 GB) DONE: 142 s compute, 399 s total
  incl. write; north-pole block = C.29 exact. Outputs cont_prop_L{1,2,4}/Dinv.0.h5.

## Key physics result
The $k$-integral + $\iota_3$ sum collapse Eq. (C.28) to a **pure-exponential double sum** (no
Bessel/Ci): $G_{11}=\sum_{m,n}\frac{1}{4\pi c^2}A A'^*\,\mathrm{sgn}(\tau)e^{-\lambda|\tau|}$ etc.,
$\lambda=n+|m|+\tfrac12$. The missing $\iota_3$ in the signed eigenvalue $\Lambda_K=\iota_3\lambda$
is what makes the diagonal odd in $\tau$ (matching the $\mathrm{sign}(t)$ of C.29).

## Open questions — all resolved by user

> All earlier open questions are now closed by user direction:
> - **A (coincident handling):** moot — only $t\ne t'$ correlators (decision 6).
> - **B (deliverable level):** the all-to-all $G(x,x')$ matrix in the `Dinv.h5` schema
>   (decision 2). Correlator-level continuum prediction is a possible later add-on, not now.
> - **Spin-connection frame / spinor basis:** not a concern — gauge-invariant quantities only.
> - **Normalization:** irrelevant (decision 7).

## Ordered implementation chunks

### Chunk 1 — angular + normalization primitives, validated standalone
Files: `src/both_3d/cont_prop_eigbasis_claude.cpp` (new)
- `xi(|m|, n, z)` via Jacobi $P^{(\alpha,\beta)}_j$. Use the three-term recurrence (stable;
  note App. C footnote: $\alpha,\beta$ here violate the usual $>-1$, but $j,\alpha+j,\beta+j$
  are nonneg integers so Rodrigues/recurrence still hold). Cross-check vs `gsl_sf_hyperg_2F1`.
- `c2(|m|, n)` via `gsl_sf_lngamma` (log-space, factorials of half-integers).
- Unit tests: edge values Eq. C.15; orthogonality $\int dz\,\xi\xi$ vs $c^2$ by quadrature.

### Chunk 2 — the $\tau$-kernels (k-integral) + temporal validation
Files: same `.cpp`
- Implement the three $\tau$-kernels by `gsl_integration_qagi`; assemble the $2\times2$
  spinor structure with the $\iota_3$ sum.
- Specialize to $\theta=\theta'=0$ and reproduce Eq. C.29 (print the $\Delta=1$ tower,
  compare to closed form). Gate: agreement to quadrature tol.

### Chunk 3 — general $G(\theta,\theta',\Delta\phi,\tau)$ with truncated double sum
Files: same `.cpp`
- Double sum over $n\in[0,40]$ and half-integer $m$ (fixed $n_{\max}=40$; $\tau\ne0$ makes it
  uniformly convergent via $e^{-M|\tau|}$). Also validate against the pure-$S^2$ CFT anchor
  $\sigma_1/(4\pi\sin(\theta/2))$ (Eqs. `psipsi_S2`/`psipsi_S2_ev`, $(1,2)$ comp).
- Self-consistency: $\sigma_3$-symmetry/parity relations from App. C; compare $\theta=\theta'=0$
  slice to Chunk 2.

### Chunk 4 — lattice driver: read coords, build all-to-all in `Dinv.h5` schema
Files: same `.cpp` + `run_cont_prop_claude.sh` (build+run handoff, tees to log)
- Read `pts_n<L>.dat`; build $(\theta_{i_x},\phi_{i_x})$ in lattice site order.
- Separation-reduced: compute $G$ over inequivalent $(\theta_{i_x},\theta_{i_y},\Delta\phi,\tau)$,
  $\tau\ne0$ only; then scatter into the **row-major $N^2$** `Dm_inv` using the decode
  $r=N_x s+N_S i_x+a$ (decision 4). Same-time diagonal blocks ($\tau=0$) left zero.
- Write `data_<esnid>/prop_exact_L<L>/Dinv.<k>.h5` with the **exact same datasets** as
  `jj_propagator_exact_claude.cu` (`Dm_inv/{real,imag}`, `N`,`N_REFINE`,`Nt`,`n_sites`,
  `parity`,`mP`,`k`,`complete`), atomic `.tmp`+rename. (`esnid` mirrors theirs, e.g.
  `free_vmRe0.000000vmIm0.000000`; choose a distinct tag so it doesn't clobber the lattice file.)
- Cross-check a few decoded entries against their lattice `Dinv.h5` (up to overall norm).

## Notes
- No CUDA; pure CPU + GSL. Cost for $L=4$ ($N_s=162$): $\binom{163}{2}\approx1.3\times10^4$
  spatial pairs $\times N_t$ separations $\times$ (double sum to $n=40$) — minutes on CPU,
  parallelizable over pairs with OpenMP (cap 4 cores when Claude runs it).
- Output is the full $N^2$ matrix to match the schema; memory same as theirs
  (L1 144 MB / L2 1.8 GB / L4 27 GB per real+imag) — for L4 stream rows to disk if needed.
- Filenames carry `_claude`; comments use LaTeX macros, no Unicode.
