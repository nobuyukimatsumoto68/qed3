# Distillation implementation plan (perambulator storage + pipeline)

Implementation plan for **exact (deterministic) distillation** on QED3 $S^2\times\mathbb{R}$ (massless overlap,
2-component GW, $U(1)$), to compute the fermionic all-to-all for the $\sigma\sigma$ four-point / coupled GEVP
with $F^2$. **NO stochastic estimators** anywhere: dense per-timeslice eigensolve + exact overlap solves per
mode -> the perambulator is exact and bit-reproducible (test X2). Method reference: Peardon 0905.2160 (exact
distillation). (The stochastic-LapH variant, Morningstar 1104.3870, is NOT used; with the option-(D) Wilson
basis it is not even "LapH".) Context: the mixing CHECK is DONE + positive
(`sigma_f2_mixing_analysis_claude.py`); distillation is for the FULL conn / coupled GEVP. See
`sigma_sigma_f2_mixing_impl_plan_claude.md`.

## Physics / goal
Replace the all-to-all propagator by the **evec-to-evec perambulator** $\tau_{kl}^{\alpha\beta}(t',t)=
v_k^\dagger(t')\,D_{ov}^{-1}\,v_l(t)$, where $V(t)=[v_1..v_{N_v}]$ are the lowest $N_v$ eigenvectors of the
$U(1)$-covariant spatial Laplacian on the fixed $S^2$ mesh. Every diagram (loops $D_S,D'_S$; meson $C_S$;
four-point $V_S,S_S,T_S$; the internal $t\to t$ return $=\tau(t,t)$) becomes a trace over the $N_v\times N_v$
matrices. Validate at L1 ($N_v=12$ = ALL modes = EXACT local), then scale to L2+.

## Key simplifications (record these)
NO $\gamma_5$-hermiticity (3D / 2-component). The replacement is the GW/parity identity of `../../qed3_v2-6.pdf`:
- **(IV.17)** (overlap unitarity): $D_{ov}^{-1}=\mathbb{1}-(D_{ov}^\dagger)^{-1}$, i.e. $D_{ov}^{-\dagger}=\mathbb{1}-D_{ov}^{-1}$
  ($\mathbb{1}$=coordinate Kronecker delta = ULTRALOCAL). **(IV.18)** parity of the $D_{ov}$ propagator up to the
  contact term; **(IV.13)** parity of $D$. These REPLACE $\gamma_5$-herm.

1. **No backward solves -- grounded in (IV.17).** The backward perambulator is the forward one minus an equal-time
   contact: $\bar\tau_{kl}(t',t)=v_k^\dagger(t')D_{ov}^{-\dagger}v_l(t)=\delta_{t't}\delta_{kl}-\tau_{kl}(t',t)$.
   So compute/store ONLY forward $\tau$; any backward $D_{ov}^{-\dagger}$ leg = contact - $\tau$, no extra solve.
   (The trivial adjoint relation $\bar\tau=\tau^\dagger$ also holds -- use it as a reality/symmetry CHECK.)
2. **FS furnished leg COLLAPSES (IV.17).** $(1-D_{ov}^\dagger)D_{ov}^{-\dagger}=D_{ov}^{-\dagger}-\mathbb{1}=-D_{ov}^{-1}$
   -- a plain FORWARD propagator. So the FS ($\tilde S$) legs of the meson/four-point may reduce to $-D_{ov}^{-1}$;
   whether $\tau'$ is still needed (vs a sign on $\tau$) is a per-diagram bookkeeping question -> pin vs v3-4 Ch.5.
3. **L1 is exact & trivial.** $N_s=10L^2+2=12$ spatial sites at L1, so $N_v=12$ = the full space = the exact
   local all-to-all, no smearing, cheap. Ideal validation vs the stochastic $D_S/D'_S$.
4. **FS extra object (if not collapsed by #2).** $\tau'_{kl}(t',t)=v_k^\dagger(t')(1-D_{ov}^\dagger)D_{ov}^{-1}v_l(t)$
   reuses the SAME forward solve $\psi=D_{ov}^{-1}v_l$ + one $D_{ov}^\dagger$ apply (forward apply, NOT a solve).

## HDF5 format (the perambulator file) -- THE core deliverable
One file per (ensemble, config): `data_<ens>/distill_Nv<Nv>/peram.<k>.h5`.

Convention (D): SPINOR modes -> spin is FOLDED INTO the mode index, so $\tau$ has NO separate $\alpha\beta$
axis. Complex arrays stored **split real/imag** (`.../real`,`.../imag`, `float64`; robust, h5py-trivial).
Row-major (C order). Mode index $k,l=0..N_v-1$; spinor component $j=N_s\cdot x+s$ ($x$=site, $s$=spin).

| dataset | shape | dtype | meaning |
|---|---|---|---|
| `/meta/Nv` `/meta/Nt` `/meta/Ns` `/meta/Nspin` | scalar | int | dims ($N_v\le 2N_s$; L1 full $=24$) |
| `/meta/config_k` | scalar | int | gauge config index $k$ |
| `/meta/M5` `/meta/at` | scalar | double | basis kernel $M5{=}-1$; $a_t$ (chunk-2 overlap) |
| `/meta/method` `/meta/ensemble` | scalar | string | "exact_distillation_wilson_basis"; ens id |
| `/evals` | (Nt, Nv) | float64 | $D_{W,2}^\dagger D_{W,2}$ eigenvalues $\lambda_k(t)$ (ascending) |
| `/V/real` `/V/imag` | (Nt, Nv, 2$N_s$) | float64 | spinor modes $w_k(t,j)$, $j{=}N_s x{+}s$ |
| `/peram/tau/real` `/imag` | (Nt, Nt, Nv, Nv) | float64 | $\tau_{kl}(t_{snk},t_{src})=w_k^\dagger D_{ov}^{-1}w_l$ |
| `/peram/tau_gw/real` `/imag` | (Nt, Nt, Nv, Nv) | float64 | $\tau'$ (FS-furnished leg; see (IV.17) collapse) |

CHUNK 1+2 (DONE+VALIDATED, `distill_peram_claude.cu`, `run_distill_L1_claude.sh`) write `/meta`,`/evals`,`/V`
AND `/peram/{tau,tau_gw}` (windowed $(W,W,N_v,N_v)$, `tsrc0/twin` in meta). L1 validation (free + 3 real cfgs):
T1a-T1d $\sim10^{-14}$; **T2a (exact-local) $\sim7\times10^{-7}$, T2b (IV.17) $\sim10^{-6}$-$10^{-5}$, T2c $\sim8\times10^{-7}$**
(all $\lesssim$ CG tol $10^{-5}$). Free-field shells (8,4,4). ~285 s/config real (mrhs-block later). GPU1.
Index order for `tau`/`tau_gw`: `[t_snk, t_src, k(sink mode), l(src mode)]` (spin inside $k,l$).
Backward legs $\bar\tau,\bar\tau'$ are NOT stored (= contact $-\tau$ via IV.17 at contraction). Elementals
$\Phi=V^\dagger W_0 V$ are NOT stored (cheap to rebuild from `/V` and the $\ell{=}0$ area weight).

**Chunking / compression:** chunk `(1, Nt, Nv, Nv)` (one sink timeslice) + gzip level 4.
**Size** per config $=2\times N_t^2 N_v^2\cdot16$ bytes (tau + tau_gw; split real/imag $=16$ B/complex; spin
now inside $N_v$, no $\times4$):
- L1 (full $N_v{=}2N_s{=}24$): $2\times128^2\cdot24^2\cdot16\approx 0.30$ GB/config (pre-gzip).
- L2 (full $N_v{=}2N_s{=}84$): $2\times128^2\cdot84^2\cdot16\approx 3.7$ GB/config -> heavy; truncate $N_v$
  and/or use the $|\Delta t|$ cutoff.
**SOURCE WINDOW (adopted, the real saver).** We only read a plateau at modest $\Delta t$ (the $F^2$-$\sigma\sigma$
cross died by $\tau\sim10$), so do NOT fill the full $128\times128$ peram. Fix a base source $t_0$ (default 0)
and a window $W$ (`--twin`, default 32); SOLVE at every source $t\in[t_0,t_0{+}W)$ (each solve gives all sinks
for free) and STORE the $W\times W$ source$\times$sink block $\tau(t',t)$, $t,t'\in[t_0,t_0{+}W)$. This is what
the GEVP needs: loops $\tau(t,t)$ + meson return leg $\tau(t_0,\Delta t)$ + $\tau(\Delta t,t_0)$ all live in the
window. Peram shape becomes $(W,W,N_v,N_v)$ with $t_0$ in `/meta/tsrc0`. Params `--tsrc0`,`--twin`.
- Cost: $W\,N_v$ solves/config (L1 $W{=}32,N_v{=}24$: 768, vs 3072 full).
- Storage: $2 W^2 N_v^2\cdot16$ B (L1: $\approx19$ MB/config, vs 0.30 GB full).
- Single origin $t_0{=}0$ drops full translation-avg, but distillation is EXACT (no stochastic noise) so this
  only costs gauge-noise reduction; the contraction can still partial-average $t_{src}\in[t_0,t_0{+}W{-}\Delta t)$
  within the window, and extra origins ($t_0{=}32,64,96$) can be added later if the signal needs it.
- SIZE $W$ FROM THE SINGLE-MESON PLATEAU: the established L1 $\ell{=}0$ scalar $C_S$ fit window is
  $dt\in[8,24]$ ($a_t{=}0.2$, `scalar_over_axial_vs_gsq_claude.py:27`; $\ell{=}1$ same, $\ell{=}2$ [12,22]).
  The two-meson $\sigma\sigma$ decays $\sim2m_\sigma$ (faster) so its usable range sits INSIDE $[8,24]$.
  -> **$W=32$** brackets $[8,24]$ with margin => default `--twin 32`. [Older full-$N_t$ $|\Delta t|\le T_{max}$
  window is subsumed by this.]
**Storage-vs-recompute:** perambulator is the EXPENSIVE cached object (the solves); store it so many operators
/ the whole GEVP contract offline without re-solving. If L2+ storage is prohibitive, fall back to
contract-on-the-fly (compute $\tau$, contract, store only the small correlator/GEVP matrices, discard $\tau$).

## $N_v$ choice
- **Spectrum first**: on $S^2$ the free $-\Delta$ has eigenvalues $\ell(\ell{+}1)$, degeneracy $2\ell{+}1$ ->
  SHELLS with gaps; cumulative $1,4,9,16,25,\dots$. Gauge lifts exact degeneracy; shells survive approximately.
  Choose $N_v$ = complete shells at a gap (icosahedral-covariant).
- **L1**: $N_v=12$ = all (exact local). No choice.
- **L2+**: empirical $N_v$-convergence scan on the observable ($\sigma\sigma$ eff mass / mixing cross) over a
  few shells on a handful of configs; pick smallest converged $N_v$. Zero-mom $\ell{=}0$ is low-mode dominated
  ($\Phi=V^\dagger W_0 V$ near-diagonal since $W_0$=const area) -> should converge with modest $N_v$.
- Two-meson at ZERO relative momentum -> $N_v$ set by a SINGLE zero-mom meson, NOT doubled. (Scattering tower
  = extra operators, not more modes; not needed for the $F^2$ mixing.)

## Contraction (offline, from the stored objects)
Rebuild $\Phi(t)=V^\dagger(t)W_0 V(t)$ ($\ell{=}0$). Each factor = a trace over the $N_v$(+spin) indices, built
from FORWARD perambulators only (v3-4 Ch.5): $D_S(t)=\mathrm{Tr}[\Phi\tau(t,t)]$;
$D'_S(t)=\mathrm{Tr}[\Phi\tau(t,t)\Phi\tau(t,t)]$; $C_S=\mathrm{Tr}[\Phi(0)\,\tau(0,t)\,\Phi(t)\,\tau(t,0)]$ (BOTH legs
forward -- $\tau(0,t)$ AND $\tau(t,0)$ both come from the full all-$(t',t)$ forward peram; NOT $\tau^\dagger$);
$S_S,T_S$ with the internal $\tau(t,t)$. Any backward $D_{ov}^{-\dagger}$ leg -> contact $-\tau$ (IV.17); FS
$\tilde S$ legs collapse to $-D_{ov}^{-1}$ (IV.17, key-simpl #2). Then weights $\{4,2,4,4,2,1,4,1,1,1\}$ -> the
coupled $\{F^2,FS{\cdot}FS,PS{\cdot}PS,FS{\cdot}PS\}$ GEVP.

## Files
- NEW `distill_peram_claude.cu` -- build $V(t)$ (covariant $-\Delta$ eigensolve, reuse `lanczos_claude.h`) +
  perambulator solves (mass-0 `OverlapWMass`, forward only) + $\tau'$ (extra $(1-D_{ov}^\dagger)$ apply) + h5 write.
- NEW `distill_contract_claude.py` (or .cu) -- read `peram.<k>.h5`, rebuild $\Phi$, assemble diagrams/GEVP.
- REUSE `overlap_wmass_claude.h` (mass 0), `lanczos_claude.h`, `s2n_simp.h`/geometry (spatial links + edge weights
  for $-\Delta$ -- CHECK if a covariant spatial Laplacian already exists, else build it).
- Handoff `run_distill_L1_claude.sh` (GPU, user runs).

## Chunk 1 -- covariant Laplacian definition (the foundational convention)
Reuse the Wilson spatial-hop ingredients (`dirac_ext.h:140`, `dirac_simp.h:355`): neighbor list `nns[x]`,
edge weight $\kappa_{xy}=2\,\text{link\_volume}_{xy}/(\ell_{xy}\,\overline{\ell})$ (= `set_kappa`), covariant
phase $e^{i\theta_{xy}(t)}$ from `U.sp(t,\{x,y\})`, lumped site measure $\bar a_x=$`dual_areas[x]`.

Covariant **stiffness** (gauge-covariant graph Laplacian, Hermitian in the PLAIN inner product):
$$
K_{xy}(t) = -\,\kappa_{xy}\,e^{i\theta_{xy}(t)}\ (y\in\mathrm{nns}(x)),\qquad
K_{xx}(t) = \sum_{y\sim x}\kappa_{xy}.
$$
The DEC/Laplace-Beltrami eigenproblem is the GENERALIZED $K v=\lambda\,M v$, $M=\mathrm{diag}(\bar a_x)$,
whose eigenvalues approximate the continuum $-\Delta$ spectrum $\lambda\simeq \ell(\ell+1)$ (shells).

Three self-consistent conventions for what we diagonalize and store as $V=[v_1..v_{N_v}]$:
- **(A) symmetric DEC** $L=M^{-1/2}K M^{-1/2}$ (Hermitian, plain inner product). Same eigenvalues as
  $M^{-1}K$ (physical $\ell(\ell+1)$ shells), and its eigenvectors $u_k$ are **plain-orthonormal**
  $u_k^\dagger u_l=\delta_{kl}$. Store $V=[u_k]$. At $N_v{=}N_s$ $V$ is UNITARY $\Rightarrow VV^\dagger=\mathbb{1}$
  (exact completeness, T1d) and $V\tau V^\dagger=P D_{ov}^{-1}P=D_{ov}^{-1}$ (exact local, T2a) -- because the
  completeness argument only needs plain-orthonormal columns, independent of the fermion measure inside $D_{ov}$.
- (B) generalized $Kv=\lambda Mv$, $v_k$ $M$-orthonormal ($v_k^\dagger M v_l=\delta$). Physical shells, but
  $\sum_k v_kv_k^\dagger=M^{-1}\neq\mathbb{1}$, so T2a needs $M$-insertions in $\tau$ -- messier.
- (C) pure graph Laplacian $\kappa{=}1$, $M{=}1$. Trivial completeness but NON-physical shells (bad truncation).

- **(D) SPINOR timeslice normal operator $D_2^\dagger D_2$.** Restrict the Wilson-Dirac operator to one
  timeslice (drop temporal hops) $=$ the existing 2D $S^2$ operator `DiracS2Simp` $D_2(t)$ (carries the
  spatial covariant hop $\kappa\,e^{i\theta}$ AND the SPIN CONNECTION $\Omega$, `dirac_simp.h:227`). Take
  $V=[w_1..w_{N_v}]$ = low eigenvectors of the Hermitian PSD $D_2^\dagger(t)D_2(t)$ (= small right-singular
  vectors of $D_2$), plain-orthonormal on the $2N_s$ spinor space.
  $$
  \tau_{kl}(t',t)=w_k^\dagger(t')\,D_{ov}^{-1}\,w_l(t),\qquad k,l=1..N_v\le 2N_s .
  $$
  Spin is FOLDED INTO the mode index -> NO separate $\alpha\beta$ index in $\tau$ (peram shape
  $(N_t,N_t,N_v,N_v)$, not $(\ldots,2,2)$). Full rank $N_v=2N_s$ (L1: **24**, not 12); there
  $VV^\dagger=\mathbb{1}_{2N_s}$, $V\tau V^\dagger=D_{ov}^{-1}$ EXACT (same anchor, spinor version).
  Gauge covariance: $D_2\to G D_2 G^{-1}$ ($G=\mathrm{diag}\,e^{ig_x}$, spin-diagonal), $\lambda$ fixed,
  $w_k\to G w_k$ (T1c). Vertices become $2N_s\times2N_s$ matrices $\Phi=V^\dagger W V$ (spin baked in, cleaner).
  **BASIS OPERATOR = timeslice WILSON $D_{W,2}$, NOT overlap.** `DiracS2Simp` IS the 2D Wilson-Dirac; it is
  cheap + ultralocal. Overlap restricted to a slice is time-nonlocal (Zolotarev of the Wilson normal op) and
  costly -> never used for the basis. The perambulator STILL inverts the massless OVERLAP $D_{ov}^{-1}$ (the
  physics: GW, (IV.17)); only the smearing basis $V$ is Wilson low modes. Basis $\ne$ propagator by design.
  **MASSLESS ARGUMENT (decisive):** $D_{ov}^{-1}$ is near-zero-mode dominated; $D_2^\dagger D_2$'s low modes are
  deflation-aligned with exactly those IR fermion modes, so truncated $N_v$ at L2 should be FAR more efficient
  than a scalar-Laplacian basis. Sub-choice: pure 2D `DiracS2Simp` (no temporal pieces; DEFAULT) vs keeping the
  temporal-Wilson diagonal shift (minor). Shells are spinor-harmonic ($j$), not $\ell(\ell+1)$.

**RECOMMENDATION (updated) = (D)** for the MASSLESS program (near-zero-mode-dominated -> best truncation), with
(A) as the vanilla-LapH fallback if the spinor basis complicates the interpolator bookkeeping. Both share the
exact-completeness anchor (D at $N_v{=}2N_s{=}24$, A at $N_v{=}N_s{=}12$) and plain-orthonormal $V$; both
diagonalize with `Eigen::SelfAdjointEigenSolver<MatrixXcd>` per (config,timeslice) (dense, tiny). Plain inner
products used consistently (source $=w_l$/$u_l$, sink $=w_k^\dagger\psi$).

## Ordered chunks
1. **Covariant spatial $-\Delta$ + eigensolve.** Build $L=M^{-1/2}K M^{-1/2}$ (convention (A) above) from
   $U(1)$ spatial links + `kappa`/`dual_areas`; dense $N_s\times N_s$ eigensolve per (config,timeslice); dump
   `/evals`,`/V`. VALIDATE shells at L1 (12 modes) and L2 (gaps), plus T1a-T1d. Files: `distill_peram_claude.cu`.
2. **Perambulator $\tau$ (+ $\tau'$).** Forward solves $\psi=D_{ov}^{-1}(v_l\otimes e_\beta)$ over all
   $(l,\beta,t_{src})$; contract $\tau=V^\dagger\psi$; $\tau'=V^\dagger(1-D_{ov}^\dagger)\psi$. Write the h5 above.
   Files: `distill_peram_claude.cu`.
3. **Contract + VALIDATE at L1.** Rebuild $\Phi$; form $D_S,D'_S,C_S$; check they reproduce the stochastic
   loops + the mixing cross (`sigma_f2_mixing_analysis_claude.py`) at $N_v{=}12$. Files: `distill_contract_claude.py`.
4. **Four-point + coupled GEVP.** $V_S,S_S,T_S$ traces; assemble $G_4$; the $\{F^2,\dots\}$ GEVP vs CP $\Delta_-$.
5. **L2 $N_v$ scan** + production. Files: `run_distill_*_claude.sh`.

## Validation / test plan (staged; anchor on ground-truth)
Distillation is DETERMINISTIC (eigensolve + solves given the config) -> exact reproducibility, and we have
three independent ground truths: (i) L1 with $N_v{=}12$ is the EXACT local operator; (ii) the stochastic
$D_S/D'_S$ loops + the positive mixing cross already measured; (iii) the free-field analytic propagator
(`project_cont_prop`, Eq C.28). Ordered checks:

**Chunk 1 -- eigensolve / covariant $-\Delta$:**
- (T1a) Hermiticity/positivity: eigenvalues real, $\ge0$; orthonormal $V^\dagger V=\mathbb{1}_{N_v}$.
- (T1b) **Free-field shells** (U=1): eigenvalues cluster into $\ell(\ell{+}1)$ shells with icosahedral-irrep
  degeneracies $1,3,5,\dots$; at L1 the 12-site (icosahedron graph) spectrum is analytically known -> match it.
- (T1c) **Gauge invariance**: apply a random $U(1)$ gauge transform to the config -> eigenvalues UNCHANGED,
  $v_k(x)\to g(x)v_k(x)$ (covariant). Confirms the covariant Laplacian is built right.
- (T1d) L1 completeness: $N_v{=}12 \Rightarrow VV^\dagger=\mathbb{1}$ (smearing = identity = exact local).

**Chunk 2 -- perambulator:**
- (T2a) **GOLD: L1 exactness.** At $N_v{=}12$, $V(t')\,\tau(t',t)\,V^\dagger(t)$ must equal a DIRECT
  $D_{ov}^{-1}$ point/unit-source solve between those timeslices, to CG tolerance. This validates the whole
  perambulator construction against a first-principles solve.
- (T2b) **(IV.17) identity.** Do ONE direct backward solve $\bar\tau=V^\dagger D_{ov}^{-\dagger}V$ and check
  $\bar\tau_{kl}(t',t)=\delta_{t't}\delta_{kl}-\tau_{kl}(t',t)$ (contact - forward). Validates the "no backward
  solves" grounding in (IV.17)/(IV.18); also check the FS collapse $(1-D_{ov}^\dagger)D_{ov}^{-\dagger}=-D_{ov}^{-1}$.
- (T2c) $\tau'$: check $\tau'=V^\dagger(1-D_{ov}^\dagger)\psi$ matches $(1-D_{ov}^\dagger)$ applied to the
  reconstructed propagator; solve residuals within tol.

**Chunk 3 -- contraction vs known data:**
- (T3a) **GOLD: L1 distilled loops == stochastic.** Distilled (exact) $D_S(t),D'_S(t)$ at $N_v{=}12$ must sit
  within the stochastic loop error bars (`corr_sigma_loops`), and TIGHTEN them (no stochastic noise).
- (T3b) **Mixing cross reproduces.** Rebuild $\langle O_F(t)[D_S^2+D'_S]\rangle_c$ from distilled loops ->
  reproduces the positive signal / same figure as `sigma_f2_mixing_analysis_claude.py`.
- (T3c) $C_S$ vs the existing single-meson conn (FNAL) -- same object, exact-vs-stochastic agreement.
- (T3d) Reality/reflection: PS diagonal correlator real, $C(t)=C(N_t{-}t)$ (or the correct reflection).

**Chunk 4 -- four-point / GEVP:**
- (T4a) **Diagram cross-checks**: the "factorizable" diagrams from the FULL contraction must equal their
  loop-built forms -- $E=C_S^2$, $F=D'_S(0)D'_S(t)$, $J=D_S^4$. Agreement pins the contraction indices/spin.
- (T4b) Time-reversal symmetries: $D(t)=C(N_t{-}t)$, $I(t)=H(N_t{-}t)$ numerically.
- (T4c) GEVP matrix Hermitian/PSD; the $F^2$ diagonal block == the standalone glue $F^2$ GEVP.

**Chunk 5 -- $N_v$ convergence (L2):**
- (T5a) Observable ($\sigma\sigma$ eff mass / mixing cross) vs $N_v$ (shells) -> plateau; the L1 $N_v{=}12$
  exact result is the anchor.
- (T5b) At L2, full $N_v{=}N_s{=}42$ (exact local) vs a truncated (smeared) $N_v$: SAME masses, different
  overlaps -- confirms smearing changes overlap not spectrum.

**Cross-cutting:**
- (X1) **Free-field analytic** (U=1): distilled $C_S$/$D_S$ vs the continuum/lattice free result
  (`project_cont_prop`, Eq C.28) -- an ABSOLUTE check independent of the stochastic data.
- (X2) Determinism: re-run -> bit-identical (no RNG), unlike the stochastic loops.

## Chunk 3 (DONE+VALIDATED, `distill_contract_claude.py`)
Vertex matched to `mult_Ylm_real(0,0)`+`accumulate_loop_raw`: $\Phi=V^\dagger W_0 V$, $W_0=\mathrm{diag}(\bar a_x/\sqrt{4\pi})\otimes\mathbb{1}$.
T3a GOLD PASS: distilled $D_S,D_S^{1mD}$ / stochastic $=1.0000$ (5 digits) on k=1/11/21; $D'_S$ estimator match
$1.0005$ (tightened; 0.03% = estimator contamination). **VERTEX PIN RESOLVED:** physical connected
$\sigma\sigma$ is $D'_{S,\text{phys}}=\mathrm{Tr}[\Phi\tau\Phi\tau]$ (TWO area vertices); the stochastic
`jj_sigma_loops` $D'_S$ computed $\mathrm{Tr}[\Phi\tau^2]$ (bare middle projector) -- MISSING one middle
$w_{00}$. Confirmed numerically $D'_{S,\text{phys}}/D'_{S,\text{est}}=w_{00}=0.295409$ EXACT at L1 (uniform
$\bar a_x$). Harmless const for the L1 mixing CHECK, but **at L2+ the middle $\Phi\neq c\,\mathbb{1}$** -> the
exact GEVP MUST use $\mathrm{Tr}[\Phi\tau\Phi\tau]$. FS: $D'^{1mD}$ uses $\tau_{gw}$ (estimator form matches,
noisier); the physical FS four-point vertex placement to be finalized with the four-point in chunk 4.
dual_areas at L1 = $4\pi/12$ uniform (icosahedron); L2+ needs actual $\bar a_x$ (dump from driver or mesh).

## Open questions
- Physical FS four-point vertex ($\tau$ vs $\tau_{gw}$ placement, two-vertex form) vs v3-4 Ch.5 -- pin at chunk 4.
- Store full $[t_{snk},t_{src}]$ or the $|\Delta t|\le T_{max}$ window? (decide from L1 correlator range).
- $N_v$ shells for L2 (from the chunk-1 spectrum).
- Chunk/gzip vs raw; store-vs-contract-on-the-fly at L2+.
- Eigenvector phase: per-config arbitrary phase cancels in physical traces; store `/V` so offline $\Phi$ uses the SAME $V$.
