# Operator-state correspondence for the $\{1,1,1,1\}$ level

## Physics / goal

In CFT the operator-state correspondence guarantees a **local (or nonlocal-kernel) bilinear
mesonic operator** $\bar\psi\,K\,\psi$ whose radial-quantized two-point function has the
$\{1,1,1,1\}$ state as its lowest pole. We want to:

1. **Redo the free limit cleanly** on $S^2\times R$ and pin the mode content.
2. **Confirm the $\{1,1,1,1\}$ quantum numbers** ($n$-labels of the two fermion legs, $\ell$, and the
   dimension $\Delta$).
3. **Construct the kernel $K$** that creates *this state and (ideally) only this state*, and understand
   whether it is (part of) the nonlocal propagator kernel $\bar\psi\,D^{-1}\psi$ (NM's hypothesis, not yet
   convinced) or a shell/radial projector.

This is the continuum resolution of the open puzzle already flagged in
`sigma_quantum_numbers_claude.md` (lines 126-131): $O_A=\bar\psi\,\tilde\tau\,\psi$ with
$\tilde\tau=D_\text{ov}^{-1}-\tfrac12$ couples to $\{1,1,1,1\}$ on the lattice, but that coupling was an
*observation*, not derived, and the $-\tfrac12$ even part "should" carry a $2E_0$ piece that is not seen.

## Free-limit spectrum (App C of `qed3_v2-6.pdf`; sources cited below)

Single free Dirac mode on $S^2$ (Eq C.16-C.18):
$$
D\,\psi_{m,n,\iota_3} = i\,\iota_3\,\lambda_{|m|,n}\,\psi_{m,n,\iota_3},
\qquad
\lambda_{|m|,n} = n + |m| + \tfrac12 ,
$$
with $m\in\mathbb Z+\tfrac12$ (antiperiodic in $\phi$), $n\in\mathbb Z_{\ge0}$ the radial (Jacobi) node
number, $\iota_3=\pm1$ the sign of the Dirac eigenvalue. On $S^2\times R$ (Eq C.22) the radial momentum
$k$ enters as $\lambda_{k,|m|,n}=\sqrt{k^2+(n+|m|+\tfrac12)^2}$.

The **shell** $\lambda=L$ ($L=1,2,3,\dots$) has energy $E_\text{leg}=L$ and degeneracy $4L$. Its
radial-orbital sub-families are $(n,|m|)$ with $n+|m|+\tfrac12=L$:

| shell $\lambda=L$ | sub-families $(n,|m|)$ | deg |
|---|---|---|
| $1$ | $(0,\tfrac12)$ | $4$ |
| $2$ | $(0,\tfrac32),\ (1,\tfrac12)$ | $8$ |
| $3$ | $(0,\tfrac52),\ (1,\tfrac32),\ (2,\tfrac12)$ | $12$ |

Free propagator (main text Eq V.3): $G(t)=\tfrac{1}{4\pi}\sum_{N\ge0}(N+1)e^{-(N+1)|t|}$ with $N\equiv\lambda-1$
(so the "leg carries $n=1$" of NM $=$ the $N=1$, i.e. $\lambda=2$, shell). **Two index conventions collide
here** and must be kept straight:
- App C radial node $n$ (in $\lambda=n+|m|+\tfrac12$),
- propagator index $N=\lambda-1$ (in $E_\text{leg}=N+1$).

### The $\ell=0$ scalar meson tower ($\Gamma=\mathbb 1$)

Bilinear $\bar\psi\,K\,\psi$, two legs $\Rightarrow \Delta=\lambda_1+\lambda_2$. The $\ell=0$
($\Gamma=\mathbb 1$, same-$\iota_3$) diagonal-shell scalars (validated in
`sigma_quantum_numbers_claude.md`, tables at lines 99-108):

| $(\lambda_1,\lambda_2)$ | $\Delta$ | operator that isolates it | free L1 effmass |
|---|---|---|---|
| $(1,1)$ | $2$ | $Q_{\lambda1}\sim\sigma_{00}$ | $2E_0=0.378$ |
| $(2,2)$ | $4$ | $Q_{\lambda2}$ | $2E_1\approx0.53$-$0.57$ |
| $(3,3)$ | $6$ | $Q_{\lambda3}$ | $2E_2\approx0.62$-$0.69$ |

where $Q_\lambda=V^\dagger\,\Pi_\lambda\,V$ is the Dirac-shell projector.

## CHUNK 1 RESULT (2026-09-17) -- the $\lambda=2$ shell is ONE $j=3/2$ multiplet; $n$ is not independent

Numerically confirmed (`state1111_freelimit_claude.py`):
- Every $\lambda=2$ mode has $j=n+|m|=3/2$. Shell degeneracy $4\lambda=8 = 2\times(2j{+}1)=2\times4$
  ($\iota_3=\pm$). **The shell $\lambda=2$ is a single $j=3/2$ single-particle multiplet.** Hence within a
  shell $n$ is locked to $|m|$: $n=\lambda-\tfrac12-|m|$. The two "sub-families" are just the
  $|m|=\tfrac32$ ($n{=}0$) and $|m|=\tfrac12$ ($n{=}1$) members of the SAME multiplet.
- The $\ell$-content of the shell-trace $\ell=0$ scalar $\rho=\sum_{a}\psi_a^\dagger\psi_a$:

  | mode set | $\ell=0$ | $\ell=2$ | verdict |
  |---|---|---|---|
  | FULL $\lambda=2$ shell | $2.257$ | $0$ | pure $\ell=0$ |
  | $(n{=}1,|m|{=}\tfrac12)$ only | $1.128$ | $0.505$ | **$\ell=0\oplus\ell=2$ mix** |
  | $(n{=}0,|m|{=}\tfrac32)$ only | $1.128$ | $0.505$ | $\ell=0\oplus\ell=2$ mix |

**Consequence for reading (B):** projecting onto $(n{=}1,|m|{=}\tfrac12)$ alone is NOT rotationally
invariant and does NOT define a pure $\ell=0$ state -- it leaks $\ell=2$. The rotation-invariant
$\Delta=4$, $\ell=0$ state is the FULL $\lambda=2$ shell (reading A). The two meanings of "leg carries
$n=1$" must be separated:
- **V.3 propagator index** $N=\lambda-1=1$ $\Rightarrow$ full $\lambda=2$ shell (reading A) -- clean.
- **App C radial node** $n=1$ $\Rightarrow$ $|m|=\tfrac12$ members only (reading B) -- not an $\ell=0$ irrep;
  only becomes physically distinct under icosahedral breaking (ties to the m-variational study, qed3-60).

Recommendation to NM: the physical $\{1,1,1,1\}$ operator is the full $\lambda=2$ shell projector
$Q_{\lambda2}$ ($\Delta=4$, pure $\ell=0$). The $n{=}1$-vs-$n{=}0$ split is the internal $|m|$-structure of
one multiplet, resolved only when rotation is broken.

**DECISION (NM, 2026-09-17): include the whole $j=3/2$ multiplet (reading A).**

## CHUNK 2 RESULT -- the operator and its dimension

The creating operator is the bilinear with the $\lambda=2$ Dirac-shell projector as kernel:
$$
O_{\{1,1,1,1\}}(t) \;=\; \bar\psi\,\Pi_{\lambda=2}\,\psi
\;=\; \sum_{a\in\lambda=2}\bar\psi_a\,\psi_a ,
\qquad
\Pi_{\lambda=2}=\sum_{a\in\lambda=2}|a\rangle\langle a| ,
$$
i.e. the shell-trace scalar over the eight $j=3/2$ modes ($m=\pm\tfrac12,\pm\tfrac32$; $\iota_3=\pm$). In the
distillation basis $K = Q_{\lambda2}=V^\dagger\Pi_{\lambda2}V$.

**Dimension (exact, analytic).** $O$ creates a fermion-antifermion pair with both legs forced into shell
$\lambda=2$. The radial-quantized intermediate energy is the sum of single-particle energies,
$E=\lambda_1+\lambda_2=2+2=4$, and the projector admits ONLY this value, so the free two-point is a SINGLE
exponential
$$
\langle O(t)\,O(0)\rangle \propto e^{-\Delta t},\qquad \Delta = 4 ,
$$
$\ell=0$ pure (Chunk 1 table: $\ell2=0$ for the full shell). This matches the free propagator V.3
($G=\tfrac1{4\pi}\sum_N(N+1)e^{-(N+1)|t|}$; the shell-2 leg decays $e^{-2|t|}$, two legs $e^{-4|t|}$).

**Lattice cross-check (already on disk).** `sigma22_lattice_claude.py` measures $Q_{\lambda2}$ on the free
lattice $\tau$: effmass $\to 2E_1$ = the $(2,2)$ state (free L1 $\approx0.53$-$0.57$ in $a_t$ units, L2
$\approx0.68$-$0.70$; the free lattice $2E_1$ shifts with refinement, continuum $\Delta=4$).

**Conclusion.** $\{1,1,1,1\}$ = the $\Delta=4$, $\ell=0$ scalar, both fermion legs in the $\lambda=2$
($j=3/2$) shell; its clean creating operator is $\bar\psi\,\Pi_{\lambda=2}\,\psi$. NM's $\bar\psi D^{-1}\psi$
guess is the RIGHT family (nonlocal, $\sigma_3$-odd, reaches $\lambda>1$) but NOT a pure creator: $D^{-1}$
weights every shell by $1/(i\iota_3\lambda)$, so it is $\{1,1,1,1\}$ plus $1/\lambda^2$-suppressed higher
shells -- the projector $\Pi_{\lambda2}$ is the exact one. (Q2 chose projector-only, so the $D^{-1}$
mechanism is left as the noted follow-up.)

## The $\{1,1,1,1\}$ identification (to confirm with NM -- see Open Questions)

Leading reading: $\{1,1,1,1\}$ is the $\ell=0$, $\Delta=4$ scalar with **both legs in the $\lambda=2$ shell**
($N=1$ on each leg), i.e. the $(2,2)$ / $2E_1$ state. The two candidate refinements of "both legs $n=1$":

- **(A) full $\lambda=2$ shell**: the $\ell=0$ scalar built from the whole shell (both $(0,\tfrac32)$ and
  $(1,\tfrac12)$ sub-families). Created cleanly by $Q_{\lambda2}$.
- **(B) the $(n{=}1,|m|{=}\tfrac12)$ sub-family only**: a specific radial-node-$1$ state inside $\lambda=2$.

The label "$\{1,1,1,1\}$" (four unit quantum numbers) suggests reading (B): each leg $= (n{=}1,|m|{=}\tfrac12)$,
two legs $\to$ four indices all $=1$ if the leg label is $(n,2|m|)=(1,1)$. **Needs NM to confirm the exact
index convention.**

## Candidate kernels $K$ for the creating operator

All act as $\bar\psi(x)\,K(x,y)\,\psi(y)$; in the mode basis a kernel with eigenvalue $\kappa_{\lambda}$ on
shell $\lambda$ weights each $(\lambda_1,\lambda_2)$ pair by $\kappa_{\lambda_1}^*\kappa_{\lambda_2}$.

| kernel $K$ | shell weight $\kappa_\lambda$ | expected two-point | note |
|---|---|---|---|
| local $\mathbb 1$ ($\sigma_{00}$) | $1$ (all shells) | ground $2E_0$ dominates | cannot reach $\{1,1,1,1\}$ |
| shell projector $\Pi_{\lambda2}$ | $\delta_{\lambda,2}$ | pure $2E_1$ | clean but "hand-built" |
| $D_\text{ov}^{-1}$ (NM) | $1/(i\iota_3\lambda)$ | all shells, $\propto1/\lambda^2$ weight | nonlocal, physical |
| $\tilde\tau=D_\text{ov}^{-1}-\tfrac12$ ($O_A$) | $1/(i\iota_3\lambda)-\tfrac12$ | observed $\to2E_1$ | $\sigma_3$-mixed, mixes w/ $\sigma^2$ |
| radial-$n{=}1$ projector | picks $(n{=}1,|m|{=}\tfrac12)$ | pure reading-(B) state | tests (B) |

Central question: does $\bar\psi D^{-1}\psi$ (or its contact-subtracted $\tilde\tau$) actually *single out*
$\{1,1,1,1\}$, or merely fail to be orthogonal to it? The $1/\lambda$ weighting does **not** project a single
shell -- so if the lattice $\langle O_A O_A\rangle$ plateaus at $2E_1$, that must come from the $\sigma_3$
parity structure (orthogonality to $2E_0$) plus the $1/\lambda^2$ suppression of higher shells, not from a
clean projection. Deriving this is chunk 3.

## Files

Reuse (do NOT rewrite; read/import):
- `free_wavefunctions_claude.py` -- App C.1 mode wavefunctions $\psi_{m,n,\iota_3}$.
- `sigma_pair_quantum_numbers_claude.py`, `sigma_ell_validate_propagator_claude.py` -- pair classification + $\tau$ validation.
- `sigma22_lattice_claude.py` -- shell projectors $Q_\lambda$ and the $2E_0/2E_1/2E_2$ ladder.
- `two_meson_gevp_OA_free_claude.py` -- the $\{1,\sigma_{00},\sigma_{00}^2,O_A\}$ free GEVP.

New:
- `state1111_freelimit_claude.py` (NEW) -- the clean redo: build the free modes, classify $\{1,1,1,1\}$,
  confirm $\Delta$, and evaluate each candidate kernel's bilinear two-point + overlap onto the $\ell=0$
  scalar tower $\{2E_0,2E_1,2E_2\}$.

## Chunks (each gated on NM go-ahead)

1. **Spectrum + classification.** Files: `state1111_freelimit_claude.py` (+ read `free_wavefunctions`).
   Tabulate shells $\lambda=1,2,3$, sub-families, degeneracies; build the $\ell=0$ scalar states; PIN
   $\{1,1,1,1\}$ = which modes, and $\Delta$. Output a short table; confirm against the App C spectrum.
2. **Kernel construction + two-point.** Files: `state1111_freelimit_claude.py`. For each
   $K\in\{\Pi_{\lambda2},D_\text{ov}^{-1},\tilde\tau,\text{radial-}n1\}$ compute the free-lattice $\tau$
   two-point effmass and the overlap vector onto $\{2E_0,2E_1,2E_2\}$. Decide which $K$ creates
   $\{1,1,1,1\}$ purely.
3. **Resolve the $D^{-1}$ mechanism.** Analytic: eigenvalue weighting $1/(i\iota_3\lambda)$, the $\sigma_3$
   parity of $\tilde\tau$, and where the $-\tfrac12$ even-part $2E_0$ overlap goes (the
   `sigma_quantum_numbers_claude.md` open puzzle). Write the derivation in this file.

## Open questions (resolve with NM before chunk 1)

- **Q1 (index convention).** What exactly does $\{1,1,1,1\}$ enumerate -- is it $(n_1,n_2)$ leg radial nodes,
  $(\lambda_1,\lambda_2)$ shells, or a $(\Delta,\ell,\dots)$ tuple? And is the target reading (A) the full
  $\lambda=2$ shell $\ell=0$ scalar, or (B) the $(n{=}1,|m|{=}\tfrac12)$ sub-family?
- **Q2 (what "the operator" should be).** A mathematically clean **pure projector** ($\Pi_{\lambda2}$, or a
  radial-$n1$ projector), or the **physical nonlocal bilinear** with a $D^{-1}$-type kernel ($O_A$-like)?
  Recommendation: construct both and compare -- the projector pins the state, the $D^{-1}$ kernel answers
  NM's convince-me question.
- **Q3 (scope).** Free limit only for now, or also carry the winning operator to the interacting L1 cache
  (Fin's `sigma2_flavorgeom_FULL_...400cfg_d1` / L1 perams `distill_Nv24`) as a validation?

## CHUNK 3 RESULT (2026-09-17) -- Fin's exact GEVP in the FREE theory (the below-threshold test)

Ran Fin's *exact* parity-even $6\times6$ flavor$\times$geometry GEVP
(`sigma2_Peven_6x6_gevp_claude.py`, ops $\{$PP,FF$\}\times\{\sigma^2_{00},O_{2m},O_{1m}\}$, block-Hankel
$\Delta t=[0,2,4]$ reb5@4 $T0=3$) but with **free perambulators** (`ENS=free`, `data_free/distill_Nv24`,
single exact config). Free reference scales: single scalar meson $m_\sigma=2E_0=0.378$; two-meson threshold
$2m_\sigma=0.756$.

Free levels (dimensionless $a_t m$):

| level | $t=7$-$10$ | identity |
|---|---|---|
| m0, m1 (doubled) | $\approx0.58$-$0.61$ | the $(2,2)$ single meson $2E_1$, parity-doubled |
| m2 | $0.757$-$0.760$ (flat) | the two-meson $2m_\sigma=0.756$ -- **exactly at threshold** |
| m3, m4 | $\gtrsim0.98$ | excited |

**The point.** In the free theory there is NO interaction, hence NO binding, yet the GEVP still shows a level
(m0/m1 $\approx0.58$-$0.61$) sitting **below** the two-meson (m2 $=0.757$). That sub-threshold level is the
$(2,2)$ single meson, which reaches $\sigma^2$ through diagram A's same-slice ($\sigma_3$-odd) contraction --
it is orthogonal to the $\Delta=2$ ground (hence the basis' lowest state is $2E_1$, not $0.378$) exactly as the
$O_A$ note predicts. The two-meson itself sits precisely at threshold $0.756$ (m2), confirming no binding.

**Resolution of NM's paradox.** In the CONTINUUM the $(2,2)$ single meson ($\Delta=4$) is degenerate with the
two-meson threshold ($2\times\Delta_\sigma=4$) -- "$\Delta=4$ is on top of $2m$." On the LATTICE, dispersion
splits the two continuum-degenerate $\Delta=4$ states, pushing the single meson $(2,2)$ BELOW the two-meson.
So a level below threshold is NOT evidence of a bound two-meson; the free GEVP produces exactly such a level
with zero interaction. Fin's interacting $0.62$ ("below $2m_{PS}=0.644$") is therefore the interaction-shifted
$(2,2)$ single meson, not a sub-threshold two-meson bound state -- as NM concluded.

(Plot: `figs/sigma2_Peven_6x6_reb5_T03_free_claude.png`. Its dashed $0.644$ line is the *interacting*
$2m_{PS}$, not the free threshold; the free threshold is m2$=0.756$.)

## Refs (algorithm / method sources)

- `qed3_v2-6.pdf` App C, Eq (C.16)-(C.22) -- free $S^2\times R$ Dirac spectrum $\lambda=n+|m|+\tfrac12$; the
  method follows Abrikosov / the Dirac-on-$S^2$ solution (Ref. [81] therein).
- `qed3int_v3-4.pdf` Fig 1 diagram A (Eq 5.5) -- the same-timeslice contraction = the nonlocal equal-time
  kernel.
- Existing project notes: `o_a_operator_note_claude.md`, `sigma_quantum_numbers_claude.md`,
  `point_split_operator_impl_plan_claude.md`, `e0e1_continuum_impl_plan_claude.md`; `project_cont_prop`
  (reconstructed free propagator).
