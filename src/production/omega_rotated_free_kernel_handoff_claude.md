# Handoff: the $\Omega$-rotated free-kernel construction of the $(2,2)$ interpolator

Author: "Fin: Two-meson" session, for **qed3-79 [745ad8] "Explore: $\Omega$-rotated free approx"**. 2026-09-18.
Self-contained: assumes no prior context from the parent thread.

Note on notation: **$\Omega$ here means the GAUGE rotation** $\Omega(x)\in U(1)$ (the thing this note is about). Our
timeslice Wilson-Dirac operator also carries a *spin connection* (unrelated); I call that "the spin connection", never
$\Omega$.

---

## 1. Purpose in one paragraph

We want a **gauge-invariant single-meson interpolator that overlaps the $(2,2)$ state** -- the $0^{++}$, $\ell=0$,
$\Delta=4$ excitation in which both constituent fermions sit in the second (spatial angular) shell $\ell=3/2$ on the
$S^2$ timeslice. Every attempt to build it from the *interacting* lattice eigenbasis fails (Sec 3). The idea to
pursue here: build the $\ell=3/2$ **projector in the FREE limit** (where the shell is exact and analytically known),
and make it gauge-invariant by **dressing it with the gauge rotation** $\Omega$ obtained from a moderate per-timeslice
Wilson flow + Landau gauge-fix:
$$
K_{(2,2)}(x,y) \;=\; \Omega(x)\,P_{3/2}^{\text{free}}(x,y)\,\Omega^\dagger(y).
$$
This sidesteps the interacting eigenbasis entirely: the shell content comes from the free projector, and $\Omega$
carries the gauge frame so $\bar\psi K \psi$ is a legitimate gauge-invariant operator. $K$ then enters the meson
correlator (Sec 9) and the $0^{++}$ GEVP.

**SANDWICH DIRECTION (corrected 2026-09-18 after qed3-79's covariance check -- see Sec 4).** With the fixing rotation
$\Omega$ defined so the SMOOTH config is $U^\Omega=\Omega U\Omega^\dagger$, the physical interacting shell projector is
$$
K_{(2,2)} \;=\; \Omega^\dagger\,P_{3/2}^{\text{free}}\,\Omega
\qquad\Longleftrightarrow\qquad K_{(2,2)}(x,y)=e^{-i\theta(x)}\,P_{3/2}^{\text{free}}(x,y)\,e^{+i\theta(y)} .
$$
(The earlier draft wrote $\Omega P^{\text{free}}\Omega^\dagger$ -- the INVERSE; that does not overlap $(2,2)$.) The one
caveat is a naming convention: if $[5c528b]$'s `landau_fix` returns the inverse-fixing rotation, then "their $\Omega$" =
"my $\Omega^\dagger$" and their `Omega P Omega^dag` already equals this -- see Sec 4/8. Lock the orientation by the
covariance algebra (Sec 4) and validation step 4 (Sec 10), NOT by the random-gauge check.

---

## 2. Physics target: the $(2,2)$ state

- Our theory: **compact $U(1)$ QED3** (2+1D) with $N_f$ two-component overlap fermions, on a **simplicial
  $S^2\times R$** lattice. Refinement $L$ = `N_REFINE`; spatial sites per timeslice $N_s = 10L^2+2$ (L1: 12, L2: 42);
  temporal extent `Nt=128`. Gauge links are $U(1)$ phases $U_i(x)=e^{i a_i(x)}$ on the $S^2$ spatial edges.
- The scalar sector: $\sigma = \bar\psi\psi$ ($\ell=0$). The two-scalar / four-fermion channel $\sigma^2$ is the
  $0^{++}$ we study (Chester-Pufu $(\bar\psi\psi)^2$ vs $F^2$ mixing, arXiv:1603.05582).
- Around threshold $2m_{PS}$ (interacting L2: $m_{PS}=0.3527$, $2m_{PS}=0.7054$) there are **three $\Delta=4$
  states**: the **$(2,2)$ single meson** ($\approx0.65$ interacting L2), the **two-meson** ($\approx0.72$), and an
  **excited** ($\approx0.80$). The $(2,2)$ is the one this construction targets as a clean single-meson interpolator.
- Free-limit reference energies (validate against these): single-meson $2E_0=m_{PS}=0.378$ (L1), first excited
  $2E_{3/2}=(2,2)=0.556$, two-meson $2\times2E_0=0.756$.

Why it matters: a good $(2,2)$ interpolator lets the coupled GEVP separate the $(2,2)$ single meson from the
genuine two-meson (both $\Delta=4$, near-degenerate), which is needed to interpret the $\sigma^2$ / $F^2$ mixing.

---

## 3. Why the interacting-eigenbasis routes fail (settled, so you don't redo them)

The distillation basis $V(t)$ = low eigenvectors of the timeslice Wilson normal operator $D_W^\dagger D_W$
(`DiracS2Simp`, $M_5=-1$; a "smearing" basis). Attempts to define the $\ell=3/2$ shell from lattice objects:

1. **From the perambulator** $M=\tau(t,t)-\tfrac12$ (contact-subtracted equal-time propagator, anti-hermitian). Its
   $|$eig$|$ ranks cluster into shells in the FREE limit ($0.131{\times}4=\ell{1/2}$, $0.078{\times}8=\ell{3/2}$,
   $0.023{\times}12=\ell{5/2}$) but **interacting the clusters dissolve** and the bare $\ell3/2$ diagonal slides all
   the way down to $m_{PS}$ (order-1 shell mixing at $g^2=1$).
2. **From `/evals` = $|D_W|^2$**: NON-monotonic in $\ell$ even free ($\ell3/2$ at $1.318$ sits *below* $\ell1/2$ at
   $1.362$; the Wilson $r$-term $\propto\ell(\ell+1)$ spoils the magnitude ordering), and interacting shows no
   clustering at all.
3. **From the Dirac eigenvalue** $D_W\xi=\lambda\xi$: the monotonic-in-$\ell$ quantity is $\mathrm{Im}\,\lambda$
   (naive part $\propto\ell+\tfrac12$) -- but that is a **free-limit** statement; interacting, rotational symmetry is
   broken and the Wilson-Dirac spectrum scatters into a 2D cloud, so no eigenvalue orders by shell.

**Conclusion (firm):** interacting there is NO kinematic shell projector -- $\ell$ is not a good quantum number once
the gauge field breaks rotational symmetry. A fixed free-limit projector inserted directly is **gauge-variant**
($\bar\psi P^{\text{free}}\psi \to \bar\psi g^\dagger P^{\text{free}} g\,\psi \neq \bar\psi P^{\text{free}}\psi$).
The $\Omega$-dressing is precisely what restores gauge invariance while keeping the (exact, free) shell content.

(A fallback that DOES work variationally -- for cross-checking your results -- is the perambulator **shell-window
GEVP**: `sigma2_shell22_gevp_claude.py`, orthogonal $M$-eigenspace rank windows $[0,4],[4,12],[12,24]$ + a GEVP that
orthogonalizes $m_{PS}$ out. It gives the $(2,2)\approx0.65$ at interacting L2. Your $\Omega$-route should reproduce
that number.)

---

## 4. The idea: $\Omega$-rotated free projector

**(a) Which sandwich = the physical projector (covariance).** Fixing rotation $\Omega$: $U^\Omega=\Omega U\Omega^\dagger$
is the smooth (Landau) config. By covariance $D_W[U^\Omega]=\Omega D_W[U]\Omega^\dagger$, and $U^\Omega\approx$ free
$\Rightarrow D_W[U^\Omega]\approx D_W^{\text{free}}\Rightarrow D_W[U]\approx\Omega^\dagger D_W^{\text{free}}\Omega$.
Eigenvectors follow the similarity, so the interacting $\ell=3/2$ eigenprojector is
$$
P^{\text{int}}_{3/2}\;\approx\;\Omega^\dagger\,P^{\text{free}}_{3/2}\,\Omega \;=\; K_{(2,2)} .
$$
This is the one that overlaps $(2,2)$. (Not $\Omega P^{\text{free}}\Omega^\dagger$.)

**(b) It is gauge-invariant (and $\Omega P\Omega^\dagger$ is not).** Under a gauge transform $h$ of the ORIGINAL
config ($U\to U^h$, $\psi\to h\psi$), the Landau gauge is unchanged, so the fixing rotation goes $\Omega\to\Omega h^\dagger$.
Then $K=\Omega^\dagger P\Omega\to(\Omega h^\dagger)^\dagger P(\Omega h^\dagger)=h(\Omega^\dagger P\Omega)h^\dagger=hKh^\dagger$,
so $\bar\psi K\psi\to\bar\psi h^\dagger(hKh^\dagger)h\psi=\bar\psi K\psi$ -- INVARIANT. Whereas $\Omega P\Omega^\dagger\to
\Omega h^\dagger P h\Omega^\dagger\neq hKh^\dagger$ -- NOT invariant. So under the fixing-rotation convention the two
sandwiches differ even in gauge invariance, and Sec-10 step 3 IS a real discriminator (contrary to a first guess).

And because a moderate Wilson flow + Landau gauge-fix brings the *frame* close to the free one, the
gauge-rotated Dirac operator $\approx$ the free operator on the smooth sector, which **justifies** using the free
$P_{3/2}$ in the first place. So the projector construction is trivial (free / continuum), and all the interacting
content is in $\Omega$.

---

## 5. The method (moderate flow + Landau gauge-fix), from agent `main-thread [5c528b]`

**Important framing correction** (this cost them time): do NOT flow to pure gauge -- that trivializes the config and
throws away the physics. Instead: **moderate** Wilson flow (smooth the UV, a few $t_0$) THEN **Landau gauge-fix**;
$\Omega$ is the Landau gauge rotation of the *flowed* config.

1. **Per-timeslice SPATIAL Wilson flow.** Flow each $S^2$ slice's spatial links independently, temporal links frozen
   (timeslices decoupled) -> one $\Omega(x)$ per timeslice, time untouched. Moderate flow time (target $\tau\sim$ a
   few $t_0$; e.g. $\epsilon\sim0.01$-$0.02$, tune `nsteps`). Stop by TARGET FLOW TIME, not a pure-gauge threshold;
   monitor plaquette rising toward 1 and the Landau functional $1-\langle\text{linkTrace}\rangle$ falling.
2. **Landau gauge-fix = extract $\Omega$.** Minimize $F[g]=-\mathrm{Re}\sum_{x,\mu}\mathrm{Tr}[g_x U_\mu(x) g_{x+\mu}^\dagger]$
   ($\equiv\|\partial_\mu A_\mu\|^2$). $\Omega$ = the minimizer $g$.
3. **Furnish:** dress the free object with $\Omega$ via the sandwich $\Omega^\dagger(\cdot)\Omega$ (Sec 9).

---

## 6. $U(1)/S^2$ specialization -- it becomes a single Poisson solve (no optimizer)

For $U(1)$, $\Omega(x)=e^{i\theta(x)}$ is abelian, so the Landau functional linearizes. With $a_i(x)=\arg U_i(x)$
(after the flow the residual is small), imposing $\nabla\!\cdot a^\Omega=0$ with $a^\Omega_\mu=a_\mu+\theta_x-\theta_{x+\mu}$
gives the **discrete Poisson equation**
$$
\nabla^2\theta \;=\; \pm\,\nabla\!\cdot a
$$
solved in one shot in the **spherical-harmonic ($Y_{\ell m}$) Laplacian eigenbasis** ($Y_{\ell m}$ = the free continuum
shell basis, so self-consistent). No gradient descent; $\theta$ unique up to a global constant + harmonic zero modes.

**Do NOT chase the sign analytically.** The $\pm$ tracks your div / $\theta$-placement discretization ($[5c528b]$'s
backward-div gives $\nabla^2\theta=+\nabla\!\cdot a$; another convention gives $-$). Instead pin it with a **one-line
convention-free check**: build $\Omega=e^{i\theta}$, form $U^\Omega=\Omega U\Omega^\dagger$, and confirm the Landau
functional $1-\langle\mathrm{Re}\,U\rangle$ (or $\|\nabla\!\cdot a\|^2$) **DROPPED** vs $U$. If it ROSE, your $\Omega$
is the inverse-fixing rotation -> negate $\theta$. Once it drops, $\Omega$ is the fixing rotation and the sandwich
$\phi=\Omega\!\cdot\!\text{in};\ P^{\text{free}}(\phi);\ \Omega^\dagger\!\cdot\!\text{out}$ lands exactly on
$K=\Omega^\dagger P^{\text{free}}\Omega$ (the physical projector). $[5c528b]$'s `landau_fix` already returns this
fixing-rotation $\Omega$ (verified in their code: $U^\Omega=\Omega U\Omega^\dagger$, $\psi\to\Omega\psi$), so if you
reuse it there is nothing to invert; the drop-check just re-confirms your $S^2$ $\theta$-solve.

**The one real subtlety -- the $2\pi$ branch / flux.** $\theta_x-\theta_{x+\mu} = -a_\mu + 2\pi n_\mu$; the integer
jumps $n_\mu$ ARE the magnetic flux (the curl of $a$). A naive $\arg$-difference Poisson solve drops the integer
jumps and returns only the trivial (flux-free) part. So the flux/winding must be handled explicitly: subtract the
vortex/harmonic piece, OR flow far enough that $|a_\mu|<\pi$ everywhere so the branch is unambiguous.

---

## 7. Topology/flux caveat, and NM's resolution

Wilson flow **preserves** topological charge / magnetic flux -- a smooth gauge frame cannot gauge away vortices or the
$U(1)$ flux sectors on $S^2$. So a naively-computed $\Omega$-dressed projector is faithful on the smooth/bulk sector
but not where flux lives; QED3 on $S^2$ genuinely has flux sectors (monopole physics is central here), so this is not
a corner case. **NM's resolution: OVERFLOW.** In their later studies, pushing the flow further ("overflow it")
handles the flux in practice, so treat the topology caveat as **not a blocker** -- flow enough (and/or handle the
branch of Sec 6) and the $\Omega$-dressed free projector is faithful. This is a key thing to verify empirically here.

---

## 8. Code map (agent `main-thread [5c528b]`, under `/mnt/baracuda_14/dwms/` -- readable, same machine)

Their code is Grid-based SU(N) plus a **coordinate-free variable-group** testbed; the variable-group files
instantiate to $U(1)$ and ARE the closest template (the $U(1)$ instantiation reduces to the Poisson solve above).

- **Per-timeslice spatial flow:** `dwf4_qcd_claude/dwf4_flow_claude.h:62` `flow_spatial(gf, eps, nsteps)` (temporal
  links frozen); staple `dwf4_qcd_claude/dwf4_su3gauge_claude.h:131` `staple_wilson_spatial` (excludes time). Also
  `:90 flow_plane` (single 2D plane), `:226 flow_holed` (open-BC defect).
- **Landau gauge-fix / $\Omega$ (BEST for us):** `dwf4_qcd_claude/dwf4_gaugefix_claude.h:66` `landau_fix<G>(U0,
  Omega, niter, alpha, verbose, tol)` -- coordinate-free Fourier-Accelerated SD, templated on group $G$; instantiate
  $G=U(1)$ and the FFT preconditioner `invhatp2 = 1/(sum_mu 2-2cos p_mu)` (the inverse lattice Laplacian, ~lines
  78-88) IS the Poisson solve. Convergence $\theta=\|\Delta\|^2/(V N_c)<$ tol; $U(1)$ converges in ~1 step.
- **Furnishing (the $\Omega$-sandwich, closest template):** `Grid/Grid/qcd/utils/FreeWilson_claude.h:283`
  `FreeLimitPreconditionerW::operator()` -- 3 lines: `phi = Omega*in;  F(phi, y);  out = adj(Omega)*y;`. Your
  $\Omega$-dressed shell projector is the SAME sandwich: rotate in by $\Omega$, apply the free $Y_{\ell m}$ $\ell=3/2$
  projector, rotate back by $\Omega^\dagger$. Per-site variable-group version: `dwf4_gaugefix_claude.h:156`
  `GaugeCovFreeWilson` (`G::mulvec` / `G::mulvec_adj`).
- **Grid SU(N) fixer (if preferred):** `Test_wilson_frameopt_claude.cc:725`
  `FourierAcceleratedGaugeFixer<PeriodicGimplD>::SteepestDescentGaugeFix(Uflowed, xform, 0.1/16.0, 1000, 1e-12,
  1e-12, true, -1, false)`; `xform`=$\Omega$. GOTCHA: `alpha` must be scaled by $1/p^2_{\max}$ (their $1/16$ in 4D;
  in **2D use $1/8$** = 1/(max lattice $\hat p^2$), or tune down). FASD ref: Davies et al., PRD 37 (1988) 1581.

Their minimal $U(1)/S^2$ recipe: (1) per timeslice, spatial links $U_i(x)$; (2) optional short spatial flow;
(3) $a_i=\arg U_i$ (branch-careful re flux); (4) solve $\nabla^2\theta=-\nabla\!\cdot a$ in the $Y_{\ell m}$ basis;
(5) $\Omega=e^{i\theta}$; dress the free $P_{3/2}$.

---

## 9. Where our lattice/objects live, and how the kernel plugs into the correlator

Our QED3 side (all under `/mnt/barracuda22/qed3/qed3/src/production`):
- **Gauge configs:** `data_<ens>/` (or the ensemble dir) `ckpoint_lat.*`; the ensemble for L2 is
  `Nf2_gsq1.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L2_hb1.000000`. Links are $U(1)$ phases on $S^2$
  spatial edges. Geometry: `geom_hopping_claude.py` (`gh.build(GEOM, L)` -> spin connection $\omega$, $\alpha$,
  nearest-neighbours `nns`, `nsite`); site positions `pts_n<L>.dat`; `antipodal_map()` in `fs_gevp_point_claude.py`.
- **Timeslice Wilson-Dirac** = `DiracS2Simp` ($M_5=-1$), 2-component spinor on $2N_s$, carrying the spatial covariant
  hop $\kappa e^{i\theta}$ AND the spin connection. (Assembled densely in `distill_peram_claude.cu` chunk-1;
  a Python re-assembly of the FREE timeslice operator would be needed to get $P_{3/2}^{\text{free}}$ on our mesh --
  or use the analytic $S^2$ spinor harmonics.)
- **Free reference:** `data_free/distill_Nv24` (L1 complete) and `data_free/distill_Nv84` (L2 complete) -- 1
  deterministic config each. Use to build/validate $P_{3/2}^{\text{free}}$ (its $\ell3/2$ shell = the $\times8$
  cluster; free effmass of that shell -> $2E_{3/2}=0.556$).
- **Perambulator** $\tau(t',t)=V(t')^\dagger D_{ov}^{-1}V(t)$ (massless overlap inverse in the distillation basis),
  stored in `data_<ens>/distill_Nv24_v2/peram.<k>.h5` (`/V (Nt,Nv,2Ns)`, `/evals (Nt,Nv)`, `/peram/tau`).

**The correlator.** The interpolator is $O_{(2,2)}(t)=\bar\psi(t)\,K_{(2,2)}(t)\,\psi(t)$ with
$K_{(2,2)}=\Omega^\dagger P_{3/2}^{\text{free}}\Omega$ a $2N_s\times2N_s$ position-space kernel per timeslice (Sec 1/4 direction). In the
distillation framework the single-meson correlator is
$$
C(t,s)=-\,\mathrm{Tr}\!\big[\,\Phi_K(t)\,\tau(t,s)\,\Phi_K(s)\,\tau(s,t)\,\big],\qquad
\Phi_K(t)=V(t)^\dagger K_{(2,2)}(t)\,V(t)\ \ (N_v\times N_v),
$$
i.e. project the $\Omega$-dressed free kernel into the distillation basis to get the mode-space vertex $\Phi_K$, then
use the existing perambulator machinery. (Compare `sigma2_shell22_claude.py` / `sigma2_meson_muweight_gevp_claude.py`
which do exactly this with a different $\Phi$.) Add $O_{(2,2)}$ to the P+ GEVP together with the $\sigma^2$
two-meson operators (see `sigma2_combined_gevp_claude.py` for the combined-GEVP scaffolding: shell block, $\sigma^2$
block, and the $\langle$single$|\sigma^2\rangle$ cross via the blob $B=A(s,t)K A(t,s)$).

---

## 10. Validation ladder (do in this order)

1. **$\theta$-solve sanity (free):** on `data_free`, $a_i=0$ (trivial links) $\Rightarrow$ $\theta=0$, $\Omega=1$,
   and $K=P_{3/2}^{\text{free}}$; its meson effmass must hit $2E_{3/2}=0.556$ (L1). This validates $P_{3/2}^{\text{free}}$.
2. **$\Omega$-extraction on an interacting config:** flow a single L2 timeslice, solve the Poisson $\theta$-equation,
   run the Sec-6 drop-check (Landau functional must DROP; else negate $\theta$) and confirm $|a|<\pi$ (branch OK) -- overflow if not.
2b. **Trivial-flux-slice sanity ($[5c528b]$'s tip):** do step 2/3 FIRST on a slice with total $S^2$ flux $=0$ and
   $|a_i|<\pi$ (unambiguous branch). Confirm the $\Omega$-dressed projector still reproduces the free $(2,2)$ there.
   Rationale: overflowing kills the flux for the drop-check, but any residual winding that survives the Poisson solve
   leaks into the LOW harmonics -- exactly where $\ell=3/2$ lives -- so it can quietly contaminate the projector.
   This step isolates "$P^{\text{free}}$ + sandwich correct" from "flux/branch handling correct." ($[5c528b]$ offered
   direct help on the branch/flux + the $Y_{\ell m}$-Laplacian swap-in -- use it.)
3. **Gauge invariance check:** apply a random $U(1)$ gauge transform to the config; $C(t,s)$ from
   $K=\Omega^\dagger P^{\text{free}}\Omega$ must be unchanged (whereas the bare $\bar\psi P^{\text{free}}\psi$ AND the
   wrong-direction $\Omega P^{\text{free}}\Omega^\dagger$ both change) -- this is a real discriminator (Sec 4b).
4. **Interacting L2 result:** the $O_{(2,2)}$ effmass / its slot in the GEVP should reproduce the shell-window GEVP
   number $(2,2)\approx0.65$ (interacting L2, `sigma2_shell22_gevp_claude.py`). Agreement = the $\Omega$-route works.

---

## 11. Open questions / deliverables for qed3-79

- Build $P_{3/2}^{\text{free}}$ on our $S^2$ mesh (analytic spinor $Y_{\ell m}$ vs the free `DiracS2Simp` eigenmodes
  -- pick whichever is cleaner; the $\times8$ shell).
- Implement the per-timeslice spatial flow + the $U(1)$ Poisson $\theta$-solve in the $S^2$ harmonic basis (our
  Laplacian, not a flat FFT). Reuse the agent's `flow_spatial` staple structure + `landau_fix<U(1)>` logic.
- Nail the flux/branch handling (Sec 6) and test NM's "overflow" claim empirically (Sec 7).
- Project $K=\Omega^\dagger P_{3/2}^{\text{free}}\Omega$ to $\Phi_K=V^\dagger K V$, feed the meson correlator, and
  compare to the shell-window $(2,2)\approx0.65$.

Ping me (Fin: Two-meson) if you need the perambulator/GEVP plumbing details or the shell-window cross-check numbers.
