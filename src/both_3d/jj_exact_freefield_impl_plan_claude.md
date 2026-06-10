# jj_exact_freefield: exact all-to-all current correlators (free field) -- impl plan

## COMMUTATOR APPROACH ABANDONED (user, 2026-06-09)

Dropped method C (commutator K from [D_ov,Theta]) and the proposed box+zero per-t K(t) builder.  Reason:
the commutator gives the UNIT-weighted slice current $\sum_n k(n,t)$; the exact/physical projected current
is $w_{tp}[n]$-weighted (Eq. 4.32).  At L=1 $w_{tp}$ is constant (vertex-transitive icosahedron) so they
match up to one constant (box+zero windowed K(0) vs exact tp_t0: residual 1.5e-5 at window W=40, validated).
But at L2/L4 $w_{tp}[n]$ VARIES; applying per-$n$ weights needs a spatially-varying $\Theta_n$, which
reintroduces spatial-link contamination (per-site $\Theta_n$ has $\delta\theta\neq0$ on every spatial link out
of $n$ -- Frobenius).  So the commutator can't cleanly give the weighted current beyond L=1 -> abandoned.
Files jj_commut_tp_claude.cc / tmp_commut_claude.sh / corr_deter_commut_L1 left in place but no longer used.
=> tp now comes from A (exact basis loop, L=1) and B (local op); D (local + continuum G) remains the CFT
test but is BLOCKED on the continuum-vs-lattice propagator convention mismatch (below).

## PENDING CHECKS (status 2026-06-09) -- what is left to run

Disk state verified: lattice L1 prop (Dinv.0.h5 +Dov) DONE; exact K (basis) data_free_Kdense_L1/K.0.h5
CACHED (tmp1) => A needs only the cheap contraction, NOT a ~50-min rebuild; corr_deter_commut_L1 DONE (C);
corr_deter_local_L1 STALE (old Omega build_W, pre-redesign) => B must re-run; cont_prop_L{1,2,4}/Dinv.0.h5
all present (continuum G for D).

L=1 RESULTS (2026-06-09, tmp_AB_claude.sh): **A (exact) and C (commutator) tp shapes overlay EXACTLY**
(1.000,0.438,0.217,0.114,0.035,0.012 at dt=1..8) => commutator method + collapse formula + antiperiodic
conn_shift all validated against the basis-loop exact K.  disc~1e-16 everywhere.  B (local) decays
DIFFERENTLY (tp 0.546,0.152,0.054; sp likewise) -- EXPECTED: the ultralocal bare-gamma current drops the
overlap resolvent time-tails, a distinct operator; not a bug, not meant to match exact at finite a.  =>
A==C done; B confirmed as the distinct ultralocal op.  NEXT = D (local + continuum G) = the CFT test.

D(L1) RESULT (2026-06-09, tmp_D_claude.sh): local op + CONTINUUM G, ran OK (disc~0).  CFT ratio
Gs/Gt is NOT -2 and NOT constant: -12.7,-4.9,-7.6,-12.2 at dt=1..4 then blows up (Gt decays very fast).
tp/sp continuum SHAPES differ (tp 1,.16,.02; sp 1,.06,.01) -- for the CFT they should share f and differ
only by -(D-1), so a flat -2.  => L=1 ratio test FAILS (~-12, erratic).  At 12 icosahedral sites the
spatial current is severely under-resolved, consistent with L=1 coarseness; does NOT yet confirm/refute
Eq. 4.28.  The L=2,4 trend discriminates coarseness vs a wrong Gs formula.  (Constants not chased, but the
~6x ratio offset AND the shape mismatch are both present.)  Output: corr_deter_local_cont_L1.

L=1 milestone (ALL runnable now; NO 50-min K rebuild needed -- K is cached):
- [x] A  exact K + lattice P: the cached K.0.h5 is CORRECT (antiperiodic BC is baked into D_ov, so K(0)
      including its wrap time-tails is right -- NO K rebuild needed).  The BC bug is ONLY in the
      contraction's time-shift: jj_contract_deter conn_shift (:131) uses tshift (:59) = plain cyclic
      (t+dt)%Nt with no -1 on wrap.  FIX = make tshift/conn_shift ANTIPERIODIC-AWARE (-1 when the time index
      wraps past Nt), reusing the cached K(0): for the free (time-translation-invariant) field
      A(t)=S_t A0 S_-t with the signed shift is EXACT.  Then re-run the cheap contraction.  (Full all-t K
      rebuild = the ~50-min+ fallback, NOT needed.)
- [x] B  local op (BARE-GAMMA) + lattice P: RECOMPILE jj_local_deter (binary jj_local_deter_L1.o is the OLD
      Omega version) + run -> overwrite the stale corr_deter_local_L1.
- [x] C  commutator + lattice P: DONE/validated (corr_deter_commut_L1; disc=2.3e-13, Gt=decay+pedestal).
- [x] cmp A/B/C at L1: tp three-way, sp two-way (A vs B); shapes (constants not chased).
- [x] D  local op + CONTINUUM G at L1: add --prop-file to jj_local_deter; run vs cont_prop_L1; compare CFT.

L=2 milestone:
- [ ] build L2 free massless propagator (dense Dov+P, N=10752, 1.85 GB) -- Titan V.
- [ ] B, C, D at L2;  cmp B/C (tp), B sp, D vs CFT.

L=4 milestone:
- [ ] build L4 free massless propagator (N=41472, 27.5 GB) -- A100-80.
- [ ] GPU-ify jj_local_deter contraction (host-dense A=W^P.P 27.5 GB won't fit at L4).
- [ ] B, C, D at L4.
- [ ] **FINAL CFT VERIFICATION: D(L4) vs Eq. 4.28 Gs / 4.31 Gt / ratio -2.**

A + B SET UP (2026-06-09): conn_shift antiperiodic fix DONE (jj_contract_deter_claude.cu: twrap_sign + signed
conn_shift, reusing cached K(0) -- exact for free).  Run both via tmp_AB_claude.sh (set CUDA_VISIBLE_DEVICES):
recompiles jj_contract_deter (A, re-contract cached K) + jj_local_deter (B, bare-gamma), writes
corr_deter_exact_L1 / corr_deter_local_L1, and prints the 3-way tp / 2-way sp shape compare.
Run: `bash tmp_AB_claude.sh 2>&1 | tee AB_claude.log`.

Deferred / code-debt (not blocking the above unless noted):
- [x] exact-path antiperiodic shift in jj_contract_deter (DONE -- signed conn_shift; no K rebuild).
- [ ] wire jj_commut_tp into run_deter_claude.sh.
- [ ] axial currents (ylm tower + axial) for the deter pipeline.
- [ ] (sidelined infra, harmless) -DNT_CLI override + --tp-only-t in kbuild from the dropped small-Nt check.

## MASTER COMPARISON PLAN -- free L1/L2/L4, verify the CFT formulas (user, 2026-06-09)

GOAL: verify the analytic CFT current-correlator formulas (Eq. 4.31 $G_t$, Eq. 4.28 $G_s$, ratio
$G_s=-(D-1)G_t=-2G_t$) in the FREE conformal field.  Four ways to get the current-current correlator
(free, U=1, Nt=128), compared by SHAPE (normalization constants differ per method -- not chased):

| method                                  | tp | sp | L feasible | role                |
|-----------------------------------------|----|----|------------|---------------------|
| A. Exact K (basis loop) + lattice P     | y  | y  | 1 only     | code check          |
| B. Local op + lattice P                 | y  | y  | 1,2,4      | code check          |
| C. Commutator + lattice P               | y  | -- | 1,2,4      | code check (sp = Frobenius-forbidden) |
| D. Local op + CONTINUUM propagator G    | y  | y  | 1,2,4      | **CFT VERIFICATION** (L=4 finest) |

KEY (user): the free continuum theory IS the CFT; its current correlator = ultralocal (e_hat.sigma)
current contracted with the CONTINUUM free propagator $G$.  So D = local-op + continuum-$G$ on the L=4
point set is the direct discretization of the CFT prediction -- THAT verifies Eq. 4.28/4.31.  A/B/C (lattice
$P$) are code-consistency checks: they must agree and approach D as $a\to0$ (the L-trend), but do not by
themselves test the CFT.

Implementation notes:
- A: jj_kbuild_exact (K built tmp1, ~2299 s) + jj_contract_deter.  L1 only.  BC: free build is t=0-only +
  conn_shift => fine at short $dt$, needs per-t/antiperiodic for large $dt`.
- B: jj_local_deter (bare-gamma; NOT yet re-run -- earlier tmp2 used the old Omega build_W).  L1/L2 host OK;
  L4 host-dense A=W^P.P ($N^2$=27.5 GB) needs the contraction on GPU.
- C: jj_commut_tp (DONE, validated L1).  Cheap O($N^2$) from dense Dov+P; L4 loads Dov+P ~55 GB host (or
  stream).
- D: jj_local_deter reads ONLY P (Dm_inv) from a Dinv-schema file + builds the current from lattice geom +
  free U => feed the CONTINUUM propagator cont_prop_L{1,2,4}/Dinv.0.h5 (other agent; same schema, on disk)
  via a propagator-path OVERRIDE (--prop-file/--prop-dir, TO ADD).  No new physics code.
- Propagators (free massless, dense Dov+P): L1 done; L2 (N=10752, 1.85 GB) Titan V; L4 (N=41472, 27.5 GB)
  A100-80.  Continuum cont_prop_L{1,2,4} already exist.

EXECUTION ORDER (chunked):
1. Re-run B at L1 (bare-gamma, lattice P).  2. Produce A at L1 (jj_contract_deter on tmp1 K).
3. Compare A/B/C at L1 (tp), A/B (sp).  4. Add --prop-file to jj_local_deter; run D at L1; vs CFT.
5. Build L2 propagator; run B,C,D at L2.  6. Build L4 propagator (A100) + GPU-ify local contraction;
   run B,C,D at L4.  7. FINAL: D(L4) vs CFT (Gs, Gt, ratio -2).

## COMMUTATOR shortcut for the projected currents from the dense D_ov (user idea, 2026-06-09)

### Correct identity (NOT a gradient-1-form / curl statement -- my earlier framing was WRONG)

`dtheta` $=\delta\theta_{wz}=\theta_w-\theta_z$ is the per-LINK phase difference (Eq. 3.25/3.26), NOT an
exact 1-form.  The exact identity (conserved_current_impl_plan_claude.md:40 + the existing check code
saved_scripts_claude/check_conserved_current_claude.cu) is, for any diagonal $\Theta=\mathrm{diag}(\theta_x)$:
$$ [D_\text{ov},\,\Theta] = \sum_{\langle w,z\rangle} \delta\theta_{wz}\, k^{wz}, \qquad \delta\theta_{wz}=\theta_w-\theta_z, $$
with $k^{wz}$ the FULL overlap current (Eq. 3.15, resolvent terms included).  Holds for ANY $\theta$
profile, **temporal AND spatial alike** -- the link type does not matter (user).  Proof = gauge
covariance: $\delta_V S=\eta^\dagger[D_\text{ov},\Theta]\xi=\sum_{wz}\delta\theta_{wz}\,\eta^\dagger k^{wz}\xi$.
check_conserved_current_claude.cu verifies it at the Wilson/$\mathcal W$ level: $[X,\Theta]\xi=\sum_\ell
\delta\theta_\ell W^\ell\xi$, $X=D_W/\lambda_M$, via `multTheta` (diagonal per-site $\theta$ multiply) +
`build_W`, over random / localized / constant $\theta$ on spatial+temporal links.

### Consequence: dense D_ov makes EVERY projected current resolvent-free

With dense $D_\text{ov}$ (already built for the propagator), $[D_\text{ov},\mathrm{diag}(\theta)]$ is an
$O(N^2)$ commutator returning the full overlap current contracted with that $\theta$ -- **no multishift**:
$$ \big([D_\text{ov},\Theta]\big)_{ij}=(D_\text{ov})_{ij}\,(\theta_i-\theta_j). $$
`kop.apply_k` runs the resolvents only because it is matrix-free; we have $D_\text{ov}$, so we skip them.
Milliseconds vs the ~38-min column-by-column `kop` build.  Applies to `sp` exactly as to `tp` -- my
earlier "sp has curl, cannot be a commutator" caveat was WRONG (artifact of misreading $\delta\theta_{wz}$).
Reuse `multTheta` + the $\theta$-profile patterns already in check_conserved_current_claude.cu.

### tp correlator collapses to slicing + a Hadamard sum -- NO matmul (user, 2026-06-09)

For the time-step $\Theta_t=\Pi_{>t}$, split each index into past ($t_i\le t$)/future ($t_i>t$).  Block form
$D_\text{ov}=\big(\begin{smallmatrix}D_{pp}&D_{pf}\\ D_{fp}&D_{ff}\end{smallmatrix}\big)$, $\Pi_{>t}=\big(\begin{smallmatrix}0&0\\0&I\end{smallmatrix}\big)$, so the current is just the off-diagonal cut blocks:
$$ K^{tp}(t)=[D_\text{ov},\Pi_{>t}]=\begin{pmatrix}0&D_{pf}\\ -D_{fp}&0\end{pmatrix} \quad(\text{slice }D_\text{ov}\text{; no multiply}). $$
Contraction with $P=D_\text{ov}^{-1}$ (massless, $D_\text{ov}P=I$):
$$ A(t)=K^{tp}(t)\,P=D_\text{ov}\Pi_{>t}D_\text{ov}^{-1}-\Pi_{>t}\equiv G_t-\Pi_{>t},\qquad G_t:=D_\text{ov}\Pi_{>t}D_\text{ov}^{-1}. $$
$G_t$ is idempotent ($G_t^2=D\Pi^2D^{-1}=G_t$), a projector similar to $\Pi_{>t}$.
- **disc** $=\mathrm{tr}A(t)=\mathrm{tr}(G_t)-\mathrm{tr}(\Pi_{>t})=0$ EXACTLY (cyclicity) -- the current one-point fn vanishes.
- **conn** $=\mathrm{tr}[A(t_0)A(t)]$.  Using $\mathrm{tr}(G_aG_b)=\mathrm{tr}(D\Pi_a\Pi_bD^{-1})=\mathrm{tr}(\Pi_{>\max(a,b)})=\mathrm{tr}(\Pi_a\Pi_b)$, the two
  projector-projector terms cancel and
$$ \mathrm{conn}(t_0,t)=2\,\mathrm{tr}(\Pi_{>\max(t_0,t)})-S(t_0,t)-S(t,t_0),\qquad S(a,b)=\!\!\sum_{t_i>a,\;t_j>b}\!\! (D_\text{ov})_{ij}\,P_{ji}. $$
$S(a,b)$ is a masked sum over the element-wise product $M_{ij}=(D_\text{ov})_{ij}P_{ji}$ (= $D_\text{ov}\odot P^{T}$) of the
two matrices ALREADY in memory.  Precompute $M$ once ($O(N^2)$) and the $N_t\times N_t$ block-partial-sum table;
then every $\mathrm{conn}(t_0,t)$ is $O(1)$ table lookups -- "retrieve row/col blocks and subtract" (user).  Whole tp
correlator = $O(N^2)$, NO matmul, NO per-insertion solve.  $\mathrm{tr}(\Pi_{>s})$ = #sites with $t_i>s$ = (Nt-1-s)\*Nx.
(Massive $P=(D_\text{ov}+m)^{-1}$: $D_\text{ov}P=I-mP$ -> the same structure with an extra $mP$ correction; free/massless
is the validation target, do that first.)

### (Anti)periodicity in the tp commutator -- TWO places (user catch, 2026-06-09)

1. **Fermion antiperiodic BC = baked into the matrix.** $D_\text{ov}$ has the $-1$ on the temporal hop at the
   $N_t-1\to0$ wrap (`dirac_ext.h:415-416`, `sign=-1` at `s==Nt-1`); dense $D_\text{ov}$ and $P$ carry it, and
   the commutator/slicing reads those entries automatically.
2. **Gauge phase $\theta$ is PERIODIC on the time circle => a step has TWO edges.** $\theta_x=\mathbb 1[t_x>t]$
   on a circle jumps $+1$ at cut $t$ AND $-1$ at the wrap.  Topological: $\oint\delta\theta=0$, so NO single
   periodic $\theta$ isolates one cut.  Hence
   $$ [D_\text{ov},\Pi_{>t}] = J(t) - J_\text{wrap}, $$
   the cut-$t$ current minus the FIXED wrap-reference current.  (At $t=N_t-1$, $\Pi_{>t}=0$ => $K_{comm}=0$:
   the two cuts coincide and cancel.)  The `kop` single-cut build avoids this by summing only slice-$t$
   temporal links (a winding $\delta\theta$, $\oint=w_{tp}\neq0$) -- not producible by a periodic $\theta$.

   **Consequence:** conn$(t_0,t)=\mathrm{tr}(A(t_0)A(t))$ for $A=K_{comm}P$ returns, by time-translation invariance,
   $$ C(t_0-t)-C(t_0-t_\text{wrap})-C(t_\text{wrap}-t)+C(0), $$
   = the physical single-cut $\langle J(t_0)J(t)\rangle$ PLUS three wrap-reference terms.  NOT identical to the
   `kop` single-cut correlator.  HANDLING (pick w/ user): (a) difference to single slice
   $[D_\text{ov},\Pi_t]=J(t)-J(t-1)$ then reconstruct; (b) keep the wrap as a fixed origin and subtract the
   three known $C$ reference terms; (c) accept the commutator object and validate vs `kop` at small $N_t$ (this
   check exposes exactly the reference terms).  [OPEN -- await user choice.]

### Temporal locality of K measured (2026-06-09) -- wrap reference is ~1e-8 at Nt=128

Cannot extract a SINGLE cut by the commutator (topology).  BUT measured the temporal locality of $D_\text{ov}$
(= $K$'s range, since $K=\partial D_\text{ov}/\partial\theta$) from the dumped dense matrix
(data_free.../prop_deter_L1/Dinv.0.h5, N=3072, Nt=128, Nx=24).  Block Frobenius norm $\|D_\text{ov}[t,t+dt]\|$
(avg over $t$), normalized to dt=0:
  dt:    1     2     4     8     16     24     32     48     64
  Dov: 2.4e-1 8.6e-2 3.2e-2 7.7e-3 7.3e-4 8.4e-5 1.1e-5 2.0e-7 5.9e-9
  P  : 2.5e-1 1.3e-1 6.3e-2 2.4e-2 4.9e-3 1.1e-3 2.3e-4 1.1e-5 7.9e-7
$D_\text{ov}$ (hence $K$) decays ~x2/timeslice early; $10^{-3}$ by dt=12, $10^{-5}$ by dt=32, ~6e-9 at
dt=64=Nt/2.  => $K$ has effective temporal support ~10-30 timeslices.  The wrap edge of $\Pi_{>t}$ sits
~Nt/2 from cut $t$, BEYOND $K$'s support, so $[D_\text{ov},\Pi_{>t}]=J(t)-J_\text{wrap}$ isolates a single cut to
~6e-9.  In conn$=C(t_0-t)-C(t_0-t_w)-C(t_w-t)+C(0)$ the wrap terms (governed by $P$'s decay) are $\lesssim
1e-6$ for $t_0,t$ away from the fixed wrap; the residual $+C(0)$ is a $t$-independent constant (likely the
$\langle Q^2\rangle$ offset seen in the stochastic tp data).  CONCLUSION: at Nt=128 commutator-tp is a valid
single-cut current to ~1e-8 if origins avoid the temporal boundary; equivalently window $K$ to $\pm$~32.

### Normalization -- NON-ISSUE (user): a single overall constant, not chased.

### RESOLVED: commutator gives tp + longitudinal sp ONLY; transverse sp is Frobenius-forbidden (user)

The commutator only ever produces weight patterns $\delta\theta_{wz}=\theta_w-\theta_z$ = a discrete EXACT
1-form (coboundary of the 0-form $\theta$).  So $[D_\text{ov},\mathrm{diag}(\theta)]$ realizes a current
projection IFF its link-weight pattern is exact (zero circulation around every cycle).
- **tp:** leaves $\{t=\text{const}\}$ are an integrable foliation; the cut weight is $d(\mathbb 1[t>t_*])$,
  exact.  On the circle $\oint\delta\theta=0$ forces the two-cut structure -> use the midpoint step
  ($dt=N_t/2=64$, two maximally-separated cuts, cross-talk ~6e-9).  Clean single-cut to ~1e-8.
- **sp:** the $\hat e^a_{wz}$ spatial direction field on $S^2$ is generically NON-integrable -- its weight
  pattern has circulation around the triangles, so it is NOT $d\theta$ for any scalar $\theta$ (Frobenius /
  de Rham).  A scalar $\theta=f(\hat x)$ gives only the GRADIENT projection $\nabla f$ (e.g.
  $\theta=Y_{\ell m}\Rightarrow\nabla Y_{\ell m}$).
- **Physical kicker:** continuity $\partial_t\rho+\nabla\!\cdot\!\vec J=0$ ties the longitudinal (gradient)
  spatial current to the charge density -> the commutator-reachable part of `sp` is already fixed by `tp`
  (no new info).  The TRANSVERSE (divergence-free) spatial current = the independent $G_s$ content = exactly
  what Frobenius forbids -> must come from the `kop` build.

PLAN:
- `tp` (and redundant longitudinal `sp`): commutator from dense $D_\text{ov}$, full $N_t$, midpoint-step
  ($dt=64$), ~1e-8.  No `kop`.
- transverse `sp` (= new $G_s$ data): `kop` build, small $N_t` for validation (cost $\propto N_t$).

IMPLEMENTED (2026-06-09):
- **Chunk A DONE**: jj_propagator_deter_claude.cu now ALWAYS dumps dense `Dov` (was `if(free_field)`;
  jj_propagator_deter_claude.cu:408).  Disk doubles (Dm_inv+Dov): L1 0.3 / L2 3.7 / L4 55 GB.
- **tp builder DONE**: jj_commut_tp_claude.cc -- standalone HOST g++ (no CUDA, no project headers).  Reads
  `Dov`+`Dm_inv` from Dinv.<k>.h5, computes conn$(t_0,t)=2\,\mathrm{tr}(\Pi_{>\max})-S(t_0,t)-S(t,t_0)$ via the
  block-sum of $M_{ij}=(D_\text{ov})_{ij}P_{ji}$ + 2D suffix sum (S).  Writes
  data_<ESNID>/corr_deter_commut_L<L>/corr.<k>.h5: `h0/tp/Gt_avg` (translation-avg), `h0/t0_<b>/tp/Vpp`
  (per-origin), `h0/disc/tp/J` (=0 sanity = PDov-I diagonal).  disc=0 analytic; conn = C(dt) + Q^2-const +
  ~1e-8 wrap.  O(N^2), no matmul.  Validation script: tmp_commut_claude.sh (builds + runs on the EXISTING
  free L1 Dinv.0.h5 which already carries Dov; prints Gt_avg + disc).  Run:
  `bash tmp_commut_claude.sh 2>&1 | tee commut_tp_claude.log`.
- NOT YET: cross-check commutator-tp vs `kop` `tp` at small $N_t`; wire jj_commut_tp into
  run_deter_claude.sh; transverse `sp` via `kop`.

## Goal / physics

Validate the analytic CFT current-correlator formulas of the note (qed3int_v2-10.pdf Sec. 4) against
the lattice, in the **free field** ($U=1$, conformal), with **NO stochastic estimator** -- exact
all-to-all traces.  Targets:
- $G_t(t)$ Eq. (4.31): $G_t=-C_J\,e^{-\Delta t}(1-e^{-t})^{-2\Delta}$  (temporal projection).
- $G_s(t)$ Eq. (4.28): $G_s=(D-1)C_J\,e^{-\Delta t}(1-e^{-t})^{-2\Delta}$  (spatial projection).
- Tower Eq. (4.34/4.35): integer-spaced dimensions, $\ell=0$ vanishing, fixed relative coefficients.
- Cross-check the **parameter-free ratio** $G_s=-(D-1)G_t$ (the thing the stochastic sp data seemed to
  violate -- see jj_corr "physical sp $\approx 0$" puzzle).  Exact traces remove noise + the single-hit
  connected bias $\mathrm{tr}(K D^{-1} K D^{-1})$, so conn and disc are clean.

## Key idea (exact trace = complete-basis estimator)

The stochastic trace $\mathrm{tr}(M)\approx\eta^\dagger M\eta$ (avg over random $\eta$) becomes EXACT if
$\eta$ runs over a complete orthonormal basis and we SUM:
$$\mathrm{tr}(M)=\sum_{i=1}^{N} e_i^\dagger M\,e_i .$$
So the existing jj estimator, with the hit loop replaced by a basis loop $\eta=e_i$ ($i=1..N$) and the
hit-average replaced by a SUM, yields exact:
- disc $J(t)=\mathrm{tr}(K(t)D_m^{-1})=\sum_i e_i^\dagger K(t)\,\phi'_i$, $\phi'_i=D_m^{-1}e_i$.
- conn $C(t_0,t)=\mathrm{tr}[K(t_0)D_m^{-1}K(t)D_m^{-1}]=\sum_i \psi_i^\dagger\,(K(t)\phi'_i)$,
  $\psi_i=D_m^{-\dagger}K^\dagger(t_0)e_i$, $\phi'_i=D_m^{-1}e_i$.
The set $\{\phi'_i=D_m^{-1}e_i\}_{i=1..N}$ IS the all-to-all propagator (its columns).  N = Comp::N
$=$ n_sites$\times N_t\times 2$ (L=1: $12\times128\times2=3072$; L=2: $42\times128\times2=10752$;
L=4: $162\times128\times2=41472$).

This reuses the *validated* jj operator algebra (op_Dm/solve_sq, ConservedCurrent K/K^dag, the
w_tp/w_sp/W_ell projection weights) unchanged -- only the source generation + accumulation differ.

## Free-field simplification

$U=1$ => $D_m$ is config-independent.  We build the propagator ONCE.  No gauge loop, no ensemble.
Massless overlap $D_{\rm ov}$ (mass=0) is the conformal vector-current operator; $K$ is mass-independent.

## LOCAL (ultralocal) current path -- proposal (Eq. 3.28), for L=2,4 feasibility

The full overlap current $K^{wz}$ (Eq. 3.27) = ultralocal $W^{wz}p_0$ + non-local resolvents
$\sum_\ell p_\ell(X^\dagger X+q_\ell)^{-1}$ (the resolvents force the multishift => the $N\times n_\text{ins}\times N_t$
blowup).  Eq. (3.28) keeps only the leading WEAK-FIELD term:
$$ J_V^{wz} \sim p_0\,\frac{A_{xy}}{\ell_{xy}}\,\hat e^a_{wz}\,\bar\Psi\gamma^a\Psi . $$
This is **ultralocal**: on link $(w,z)$ it is the point-split $\hat e^a_{wz}\sigma^a$ between the two
endpoints -- a **site-wise Pauli multiply** (cf. `valence_claude.h::mult_sigma1/2/3`, which apply
$\sigma^{1,2,3}$ at each $(s,ix,\text{spin})$), NOT a COO/`build_W`.  Coefficients dropped (only the
$\hat e^a_{wz}\sigma^a$ structure matters for the shape/ratio test); "applicable to the temporal case
too" = use the same $\hat e\cdot\sigma$ with the temporal link's direction.

Why feasible to L=4: the local current is sparse/site-local, so $K_\text{loc}^P(t)\,P$ is applied as a
cheap site-wise scatter ($O(N)$ per current, no multishift, no $N_t$ blowup, no dense $K$).  Memory =
$P$+$A$ only (55 GB at L=4 on A100-80).

Plan: keep the EXACT paths (stochastic estimator + exact all-to-all).  ADD a LOCAL-current path:
- L=1 free: run BOTH exact-K and local-K all-to-all -> compare in the notebook (validate the proxy).
- L=2,4 free: local-K only (exact infeasible), with the EXACT propagator.
Implementation (TBD): a `K_local` apply built from the point-split $\hat e^a_{wz}\sigma^a$ (site-wise,
mult_sigma-style) replacing `op_K`/`apply_k_ms` in the K-builder; contraction unchanged downstream.

### CORRECTION (user): local current = bare e_hat.sigma, NO Omega, do NOT use build_W

The local point current is $\hat e^a\,\bar\psi\sigma^a\psi$ -- only $\hat e\cdot\sigma$ = the bare
`gamma(ix,iy)=cos(alpha)sigma1+sin(alpha)sigma2` (dirac_simp.h).  It does NOT carry the spin connection
$\Omega$ (= parallel transport between sites, part of the hopping).  So the local current must apply
`gamma` DIRECTLY, NOT via `build_W` (which dresses with `Omega` + gauge link).

### DONE (2026-06-09): bare-gamma redesign of jj_local_deter_claude.cu

Replaced the `kop.build_W` projected-current build with file-scope bare-gamma builders:
- `local_W_sp(en,DW,U,t,lk)` -- spatial link: `0.5 kappa gamma(ix,iy) i exp(i u_sp)` (pos at ix) and the
  neg-oriented term at iy.  = `DiracExt::d_coo_format(pair<int,BaseLink>)` with **$\Omega\to 1$** and the
  Wilson **$-r\sigma_0$ dropped** (vector current = the gamma part only).
- `local_W_tp(en,DW,U,t,ix)` -- temporal link: the `sigma_3` part of the temporal hopping
  (`= d_coo_format(pair<int,Idx>)` with $-r\sigma_0$ dropped; no $\Omega$ in the temporal hop anyway),
  with the `%N` wrap carrying the t->t+1 hop across the time boundary.
- `build_Wp_local(proj,t,DW,U,base,w_tp,w_sp,Wp)` -- sums w_tp over dual sites (tp) / w_sp over links (sp).
Accessors: `DW.bd.gamma`, `DW.bd.kappa`, `DW.kappa_t`, `DW.sigma[3]`, `DW.lattice.map2il`, `U.sp/U.tp`.
Removed the overlap `Fermion D`, `ConservedCurrent kop`, `D.update(U)` (no longer needed -- pure Wilson
hopping).  All lambdas de-lambda'd to file-scope `static` functions (apply_Wp/trace/tr_AB/load_mat/
write_corr/write_vec) per the no-lambda style rule.  Build per-t for BOTH free and interacting (the
free cyclic time-shift breaks at the antiperiodic time boundary; bare-gamma build is cheap).  NOT yet
compiled/run -> user runs tmp2.  AXIAL ((1-D_ov)->1) = next.

### DECISION (2026-06-09): exact-path BC fix = antiperiodic-aware shift (not per-t)

The exact path (`jj_contract_deter` `conn_shift` + `jj_kbuild_exact` free=t=0-only) has the SAME bug: the
cyclic `tshift` ignores the antiperiodic $-1$ on time wrap.  FIX for exact = **antiperiodic-aware shift**
(reuse the single t=0 build, $A(t)=S_t A(0) S_{-t}$ with a $-1$ on any spinor amplitude whose time index
wraps past $N_t$).  In the FREE case K and P are bulk-translation-covariant under antiperiodic BC, so this
is EXACT (not an approximation) and cheap -- vs the 128x cost of a per-t rebuild of the EXACT overlap K
(`n_ins x N` overlap applies each).  Interacting exact would still need per-t.  `jj_local` stays per-t
(its build is cheap).  TODO: implement after the local redesign is validated; recompile the exact
propagator binary (it was built pre-stream-fix).

### RESOLVED (2026-06-09): jj_propagator_deter LU race -- FIXED (missing cudaStreamSynchronize)

The broken P (off-diag ||P[t,t+dt]||=0 for t<~26, not translation invariant) was a RACE: cuSOLVER
getrf/getrs run on a NON-BLOCKING `stream`, but the D2H copies of d_info/d_B were on the DEFAULT
stream, which does NOT wait -> the inverse was read mid-solve (block-structured garbage).  Diagnosed by
dumping dense D_ov (translation-invariant => build OK) and checking D_ov.D_ov^{-1} != I (=> LU/copy
bug).  FIX: `cudaStreamSynchronize(stream)` before each D2H read in invert_shift.  Now self-checks
in-program: `[check] || D_m . D_m^{-1} - I ||_F = 9.4e-14`, and Dm_inv is translation invariant
(||Dm_inv[t,t+1]||=0.621 const, matching numpy + the continuum reference cont_prop_L1).  The conn=0 was
purely this broken P (NOT the local code / time-shift).  NEXT: rerun the local current (now reads the
correct P) -- conn should decay.

### Local-current construction (earlier notes; superseded re Omega by the CORRECTION above)

- The $\hat e^a_{wz}\sigma^a$ spin structure = the **gamma-part of the Wilson hopping** (dirac_simp.h):
  `gamma(ix,iy) = cos(alpha)*sigma[1] + sin(alpha)*sigma[2]` (the spatial link's $\hat e\cdot\sigma$ via the
  link angle `alpha`), dressed by `Omega(ix,iy)` (spin connection) and the gauge link `exp(i u)`.  The
  Wilson `-r*sigma[0]` term is NOT the current (vector current = the gamma part only).  Temporal link:
  the analogous time-direction gamma.
- Applied **site-wise** as in `valence_claude.h::mult_sigma1/2/3` -- and this is exactly the pattern
  `meson_pq_wall_v2_claude.cu` uses: `src.mult_sigma(ab)` on the source/sink, then the spatial linear
  combination (the `gamma` weights).  So K_local reuses `mult_sigma` + the `gamma`/`Omega` link data;
  NO COO, NO overlap, NO multishift.
- **Axial local current:** the GW factor $(1-D_\text{ov})\to 1$ because $D_\text{ov}=O(a)$ (dropped in
  the leading approximation).  So the axial local current is the same site-wise sigma structure with NO
  overlap dressing -> feasible just like the vector.
- Coefficients dropped throughout (only the spin/shape structure matters for the 4.28/4.31/ratio test).
- => fully ultralocal: $K_\text{loc}^P(t)P$ via site-wise scatter, $O(N)$ per current; feasible to L=4.

IMPLEMENTATION (DONE, vector tp/sp; not yet compiled): `jj_local_deter_claude.cu` -- loads P, builds the
SPARSE projected current W^P(t) from `kop.build_W` (= the d_coo_format hopping derivative; exact
sigma/Omega/temporal, NO overlap/multishift), applies it to P by sparse row-scatter (A[i,:]+=v P[j,:]),
then disc=tr(A) / conn=tr(A0 . S_dt A0 S_-dt) free, tr(A(t0)A(t)) interacting.  Output
`data_<ESNID>/corr_deter_local_L<L>/`.  Run via `run_deter_claude.sh --local` (compiles propagator + local,
runs propagator -> jj_local per mass; skips exact kbuild/contract).  AXIAL ((1-D_ov)->1) = next.
build_W uses a scoped `COO<N>` (default ctor + dtor confirmed -> no device leak); `D.update(U)` sets
lambda_max for build_W (no apply_k).  L=4: conn_shift O(N^2 Nt) host loop should move to GPU (note).

## RESOLVED decisions (user)

- **Method**: materialize the dense $D_{\rm ov}/D_m^{-1}$ and SAVE to disk; contraction reads it back.
- **L (N_REFINE)**: compile-time (N is templated throughout); pass via a `-D` macro and a run `.sh`
  that builds one binary per L (1,2,4).  No runtime L.
- **Masses**: vector (massless $D_{\rm ov}$) AND the axial $(1-D_{\rm ov})$ GW tower.
- **Combine the propagator build with reweighting_R** (user): both already build the dense $D_{\rm ov}$
  ($N$ applies, the dominant cost).  R runs cuSOLVER `geev` (eigenvalues -> Eq. 2.5); the propagator
  runs `getrf`/`getri` (LU inverse) on the SAME matrix.  So ONE dense build per config emits both:
    - eigenvalues + R (parity only) -- the existing R output, and
    - $D_m^{-1}$ (+ $\tilde D^{-1}$ for parity) -- the all-to-all propagator for the contraction stage.
  Memory: dense $D_{\rm ov}$ + workspace $\approx 2$-$3\times N^2$ (L1 ~0.4 GB, L2 ~5 GB, L4 ~80 GB ->
  L4 won't fit one GPU; L1/L2 fine).
- **Split**: program (1) = the combined spectrum+propagator builder (extend R or a new sibling -- TBD);
  program (2) = `jj_contract_deter_claude.cu` loads the propagator and contracts into the correlators.
- **Propagator program = ONE mass per run** (CLI `--mass-re/--mass-im`).  Build the dense $D_{\rm ov}$
  ONCE (the expensive $N$ applies), then form + LU-invert the operator(s) needed for that mass and save:
  - massless ($m=0$): save $D_{\rm ov}^{-1}$ only (start here).
  - $m_F$ (real $m$): save $D_m^{-1}=(D_{\rm ov}+m)^{-1}$.
  - $m_P$ (imag $m$): save BOTH $D_m^{-1}$ and $\tilde D^{-1}=(D_{\rm ov}+m/(1-m))^{-1}$.
  Since $D_m=D_{\rm ov}+c\,\mathbb 1$, every mass is a cheap shift+LU on the one dense $D_{\rm ov}$.
  Output H5: `data_<esnid>/prop_deter_L<L>/Dinv.h5` (or similar), atomic write (.tmp+rename).
- **H5 fix first** (atomic write + recover corrupt file), then the exact work.

## L=4 free-field on A100-80 (cluster) -- feasibility + required trims

Cluster GPU = A100-80GB.  L=4: N=41472, dense $N\times N$ complex = 27.5 GB.
- **GPU**: LU needs $A$ + $B$ (RHS->inverse) = 2 x 27.5 = 55 GB; getrf panel workspace + ipiv are small.
  Fits in 80 GB.  (Skip `--with-R` at L=4: geev's $d\_VL,d\_VR$ would add 2 x 27.5 GB and blow the budget;
  R is parity-only anyway, and the free conformal test is massless.)
- **Free case = SINGLE config** (U=1): one dense build (N=41472 applies) + ONE LU (massless => one
  inverse $D_{\rm ov}^{-1}$).  No ensemble loop, no resume needed.
- **Disk**: `Dinv.0.h5` (massless) = $16N^2$ = 27.5 GB.
- **HOST-MEMORY TRIM (REQUIRED for L=4; do before running L=4 -- chunk C1b):** the current C1 host
  buffers are too fat at L=4 -- `A0` (27.5) + `B_h` (27.5) + `I_h` identity (27.5) + output `re`/`im`
  doubles (2 x 13.7) ~ 110 GB host.  Trim: build the identity ON DEVICE (memset + a tiny diagonal
  kernel, drop the 27.5 GB `I_h`); copy the device inverse straight into the row-major `re`/`im`
  output (drop `B_h`).  That brings host use to ~55 GB (A0 + output), fine on a cluster node.
  L=1/L=2 are unaffected by the current (fat) buffers, so C1b is an L=4 enabler, not a correctness fix.
- **Contraction (C2) at L=4**: loads the 27.5 GB propagator (GPU or host); K-applies are sparse/local,
  cheap.  Fits on A100-80.

## Cost comparison (L=1, N=3072), overlap (multishift) applies

- Exact dense build of $D_{\rm ov}$ = $N=3072$ applies (one per column) + $O(N^3)$ cuSOLVER LU invert
  (negligible at this N).  Zero variance (= infinite stochastic statistics).
- Stochastic Z2 (current jj), per hit $\approx 176$ outer CG solves $\times\,n^\text{out}_\text{CG}$
  iterations $\approx$ few$\times10^3$ applies; disc needs many hits to converge.
- => Exact $\approx$ 1-2 hits in apply-count, but noise-free.  Scales as $N$ applies + $N^3$ LU + $N^2$
  storage: fine for L=1,2; heavy for L=4 (41472, ~27 GB).  CONFIRM with a timing once built.

## Open questions (resolve before coding)

1. **Exact-trace method**:
   (a) [recommended] complete-basis loop reusing jj machinery (N solves; the propagator columns are
       streamed, never materialized as one big matrix) -- minimal new code, no large allocation; OR
   (b) materialize the dense $D_m^{-1}$ ($N\times N$) and do matrix contractions / save to disk.
       L=1: 144 MB; L=2: 1.8 GB; L=4: 27 GB (the heavy one).
2. **Scope of L**: L=1 and L=2 first (3072 / 10752 solves), add L=4 (41472 solves, ~27 GB if
   materialized) after the formula check passes?  Or all three now.
3. **Masses**: massless $D_{\rm ov}$ only (the conformal vector test, 4.28/4.31)?  Or also include the
   axial $(1-D_{\rm ov})$ GW legs / m_F / m_P?  (massless vector is enough to settle the 4.28 question.)
4. **Save propagator**: write $\{\phi'_i\}$ (or dense $D_m^{-1}$) to H5 so the correlator stage can be
   re-run without re-solving?  Only worth it for L=2,4.
5. **Output**: write correlators in the EXISTING jj h5 layout (keys h0/t0_b/{tp,sp,ylm,...}) so the
   SAME notebook (jj_corr_massless) reads them with one config and exact (zero-noise) values?

## NAMING CONVENTION (user, 2026-06-09) -- two orthogonal axes; "exact" is RESERVED for the Kop

- Kop axis:        `exact` (full overlap conserved current K) | `local` (ultralocal e_hat.sigma).
- inversion axis:  `stoch` (stochastic estimator) | `deter` (deterministic all-to-all D^{-1}).
So the existing jj_corr = (stoch, exact).  The new all-to-all pipeline = `deter`, and within it the
current is `exact` or `local`.  RENAME (pending, after the readability cleanup):
- jj_propagator_deter -> jj_propagator_deter   (deterministic inversion; current-agnostic)
- jj_contract_deter   -> jj_contract_deter     (deterministic contraction; loads exact OR local K)
- jj_kbuild_exact     -> jj_kbuild_exact        (KEEP: builds the EXACT K) -- pairs with jj_kbuild_local
- jj_local_deter      -> the LOCAL Kop path (local current); fold into the deter contraction
  (build local K + apply); name e.g. jj_kbuild_local / jj_contract_deter --local.
Run dir tags likewise: corr_deter_exact / corr_deter_local (was corr_exact / corr_local).

## ARCHITECTURE (implemented): 3 programs + orchestrator (user-directed)

K is mass-INDEPENDENT -> its own program + its own mass-free dir, built once per (config,L), reused
across masses.  All overlap/GPU work is in stages 1a/1b; stage 2 is pure cuBLAS + HDF5.
- `jj_propagator_deter_claude.cu` (1a) -> `data_<ESNID>/prop_deter_L<L>/Dinv.<k>.h5`  (D_m^{-1}, LU, mass-dep)
- `jj_kbuild_exact_claude.cu`     (1b) -> `data_<TAG>_Kdense_L<L>/K.<k>.h5`            (K^P(t) dense, MASS-FREE;
      free => t=0 only [contraction time-shifts]; interacting => all t.  datasets <proj>_t<t>)
- `jj_contract_deter_claude.cu`   (2)  -> `data_<ESNID>/corr_deter_exact_L<L>/corr.<k>.h5`   (A=K.P matmul; disc=tr(A);
      conn=tr(A0 . S_dt A0 S_-dt) free / tr(A(t0)A(t)) interacting.  NO overlap/geometry includes.)
- `run_deter_claude.sh` orchestrates: compile-per-L (-DN_REFINE_CLI) + 1b(once) -> 1a(per mass) -> 2(per mass).
- Layout (dirac_ext.h:65): idx = Nx*t + (NS*site+spin), Nx=NS*N_SITES => time OUTER; S_dt = rotate Nx-blocks.
- STATUS: all written; vector tp/sp done; ylm tower + axial = next.  Not yet compiled by Claude.

## Files (new; none of the production code touched)

- `jj_propagator_deter_claude.cu` -- STAGE 1: dense $D_{\rm ov}$ build + LU-invert the mass's
  operator(s) + atomic save `data_<esnid>/prop_deter_L<L>/Dinv.<k>.h5`.  Sibling of reweighting_R
  (copied; R untouched).  `--with-R` reuses the matrix for Eq. 2.5.  `N_REFINE_CLI` macro for L.
- `jj_contract_deter_claude.cu` -- STAGE 2: load `Dinv`, apply K (ConservedCurrent), contract exact
  tp/sp/ylm (+ axial $(1-D_{\rm ov})$ tower); write jj-format h5 for the notebook.
- `run_prop_exact_claude.sh` -- builds one binary per L (`-DN_REFINE_CLI=1/2/4`) and runs stage1+2.
- reuse headers: valence/overlap_wmass/conserved_current/matpoly/blocked_mat + s2/dirac stack.

## Ordered chunks

- **C1 -- propagator builder [DONE, not yet compiled]**: `jj_propagator_deter_claude.cu` -- dense
  $D_{\rm ov}$ (N applies via Op.from_cpu), cuSOLVER `Zgetrf`/`Zgetrs` LU-invert $D_m$ (+ $\tilde D$ if
  parity), row-major save, atomic write; `--with-R` geev path.  `tmp_claude.sh` builds+smoke-tests
  L=1 free massless.
- **C1b -- L=4 host-memory trim [TODO before L=4]**: device-side identity + drop `B_h`/`I_h`
  double-buffers (see L=4 section).  Needed only to run L=4; L=1/2 unaffected.
- **C2 -- contraction (exact correlators)**: `jj_contract_deter_claude.cu` -- load $D_m^{-1}$, apply
  K(a,t), accumulate disc $J_{\rm tp/sp/ylm}(t)$ and conn $C(t)$ with w_tp/w_sp/W_ell weights; physical
  $-C_{\rm conn}+C_{\rm disc}$; axial $(1-D_{\rm ov})$ legs.  Write jj-format h5 (single config).
- **C3 -- run script + compare**: `run_prop_exact_claude.sh` per-L; notebook/script overlays $G_t,G_s$
  vs analytic 4.28/4.31 and the parameter-free ratio $G_s/G_t=-(D-1)$; ylm tower vs 4.34/4.35.

## Separate small item -- H5 corruption on interrupted write

User hit "error opening a .h5" after quitting mid-write.  Prevent via **atomic write**: write each
file to `<name>.h5.tmp`, then `std::filesystem::rename` to `<name>.h5` only AFTER the 'complete'
sentinel + close.  rename(2) is atomic on POSIX, so a reader never sees a partial file and an
interrupted run leaves only a stale `.tmp` (ignored).  Apply to jj_corr_block_t / jj_disc_postproc /
reweighting_R writers.  Also: recover the current corrupt file (identify + delete so resume recomputes).
This is independent of the exact-propagator work; can do first.
