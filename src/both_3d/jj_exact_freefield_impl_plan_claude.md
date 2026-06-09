# jj_exact_freefield: exact all-to-all current correlators (free field) -- impl plan

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
  program (2) = `jj_contract_exact_claude.cu` loads the propagator and contracts into the correlators.
- **Propagator program = ONE mass per run** (CLI `--mass-re/--mass-im`).  Build the dense $D_{\rm ov}$
  ONCE (the expensive $N$ applies), then form + LU-invert the operator(s) needed for that mass and save:
  - massless ($m=0$): save $D_{\rm ov}^{-1}$ only (start here).
  - $m_F$ (real $m$): save $D_m^{-1}=(D_{\rm ov}+m)^{-1}$.
  - $m_P$ (imag $m$): save BOTH $D_m^{-1}$ and $\tilde D^{-1}=(D_{\rm ov}+m/(1-m))^{-1}$.
  Since $D_m=D_{\rm ov}+c\,\mathbb 1$, every mass is a cheap shift+LU on the one dense $D_{\rm ov}$.
  Output H5: `data_<esnid>/prop_exact_L<L>/Dinv.h5` (or similar), atomic write (.tmp+rename).
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

## Files (new; none of the production code touched)

- `jj_propagator_exact_claude.cu` -- STAGE 1: dense $D_{\rm ov}$ build + LU-invert the mass's
  operator(s) + atomic save `data_<esnid>/prop_exact_L<L>/Dinv.<k>.h5`.  Sibling of reweighting_R
  (copied; R untouched).  `--with-R` reuses the matrix for Eq. 2.5.  `N_REFINE_CLI` macro for L.
- `jj_contract_exact_claude.cu` -- STAGE 2: load `Dinv`, apply K (ConservedCurrent), contract exact
  tp/sp/ylm (+ axial $(1-D_{\rm ov})$ tower); write jj-format h5 for the notebook.
- `run_prop_exact_claude.sh` -- builds one binary per L (`-DN_REFINE_CLI=1/2/4`) and runs stage1+2.
- reuse headers: valence/overlap_wmass/conserved_current/matpoly/blocked_mat + s2/dirac stack.

## Ordered chunks

- **C1 -- propagator builder [DONE, not yet compiled]**: `jj_propagator_exact_claude.cu` -- dense
  $D_{\rm ov}$ (N applies via Op.from_cpu), cuSOLVER `Zgetrf`/`Zgetrs` LU-invert $D_m$ (+ $\tilde D$ if
  parity), row-major save, atomic write; `--with-R` geev path.  `tmp_claude.sh` builds+smoke-tests
  L=1 free massless.
- **C1b -- L=4 host-memory trim [TODO before L=4]**: device-side identity + drop `B_h`/`I_h`
  double-buffers (see L=4 section).  Needed only to run L=4; L=1/2 unaffected.
- **C2 -- contraction (exact correlators)**: `jj_contract_exact_claude.cu` -- load $D_m^{-1}$, apply
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
