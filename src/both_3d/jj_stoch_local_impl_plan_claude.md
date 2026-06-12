# Stochastic LOCAL current correlator -- implementation plan (handoff 2026-06-11, resume 2026-06-12)

## Goal
Implement the **stochastic** estimator of the LOCAL (ultralocal, bare on-site $\sigma^a$) vector-current
correlator, in the FREE and interacting theory. Scope set by the user:
- **loc only** (the on-site $\sigma$ current). **NO disp.** (K already exists in `jj_corr_block_t`.)
- channels **sp, tp, ylm**.
- **both conn and disc**, but **a SEPARATE FILE for each** (one conn `.cu`, one disc `.cu`).
- **kept SEPARATE from the K codes** -- do not fold into `jj_corr_block_t_claude.cu`.
- modeled on / "following" **`meson_pq_wall_v2_claude.cu`** (style + stochastic machinery).

## References (read these first when resuming)
- **`meson_pq_wall_v2_claude.cu`** -- the STYLE template (stochastic meson 2pt w/ wall source). Loop is
  `for t0: for h: eta.fill_z2_wall_source(rng,t0); src2=eta; op_DHD.solve -> Dinv_src2;
   for ell,em,ab: src1=eta; src1.mult_sigma(ab); op_DHD.solve -> Dinv_src1;
   Dinv_src2copy.mult_Ylm_real(ell,em,base).mult_sigma(ab);
   C(s)=sum_s psi1.slice^dag psi2.slice; write h5 ell/em/ab/h/t0`. (lines ~449-487.)
- **`jj_local_deter_claude.cu`** -- the loc CURRENT definition (deterministic, dense $P=D_m^{-1}$):
  - insertion `local_W_sigma`: bare Pauli $\sigma^a$ ($a=1,2,3$) at dual site $n$, timeslice $t$; a pure
    SITE object (no link hop, no $\kappa$, no $i$, no gauge phase, no $\Omega$, no Wilson $-r\sigma_0$).
  - channels `s1=<\sigma_1\sigma_1>`, `s2=<\sigma_2\sigma_2>` (spatial, $G_s=s_1+s_2$, Eq.4.28),
    `s3=<\sigma_3\sigma_3>` (temporal $G_t$, Eq.4.31).
  - conn (INSERTION-DIAGONAL, **same site $n$ at both times**, summed over $n$):
    `conn(t0,t) = sum_n w_site[n] tr[ \sigma^a(n,t0) P \sigma^a(n,t) P ]`.
  - disc: `disc(t) = sum_n w_site[n] tr[ \sigma^a(n,t) P ]`.
  - weight `w_site[n] = base.dual_areas[n]` -- **NO $1/\kappa^2$** (bare $\sigma$ carries no $\kappa$);
    `write_corr` folds `1/(4\pi)`.
  - ylm tower (Eq.4.35/4.36): `Sigma_lm(t)=sum_n A_n Y_lm(n^) \sigma_3(n,t)`, `A_n=dual_areas`;
    `g_l=(1/(2l+1)) sum_m tr[Sigma_lm(t0) P Sigma_lm(t) P]`, `build_Sigma_ylm` helper.
  - analysis-side (`comp_trio`): loc `G_s=(s1+s2)/(D-1)`, `G_t=s3`; ratio $\to -1$ (loc is t-indep $-1$).
- **`jj_corr_block_t_claude.cu`** -- the K stochastic conn+disc, for the DIAGONAL stochastic trick:
  - `eta.fill_z2_source(rng)` (FULL-VOLUME z2, not a wall), per hit (line 438).
  - `phi' = D_m^{-1} eta` (one solve, shared) (line 461-462).
  - conn source leg `psi_a(t0) = D_m^{-dag} K^dag(a,t0) eta`; sink `K(a,t) phi'`;
    `conn(a,t) = psi_a(t0)^dag [K(a,t) phi'] = eta^dag K(a,t0) D^{-1} K(a,t) D^{-1} eta`
    `=> E[conn(a,t)] = tr[K(a,t0) D^{-1} K(a,t) D^{-1}]` -- **insertion-DIAGONAL, exact per $a$**
    (same eta both ends isolates the same insertion). One source solve PER insertion $a$ per $t_0$.
  - disc trace `J(t) = sum_a w_a eta^dag K(a,t) phi'` rides the same sink pass (RAW; 1/4pi at analysis).
  - weights `w_tp[n]=dual_areas/kappa_t^2`, `w_sp[il]=link_volume/kappa^2` (K carries $\kappa$; loc does NOT).
- helpers in `includes/valence_claude.h`: `fill_z2_source`, `fill_z2_wall_source(rng,s)`, `slice(s)`,
  `mult_sigma(a)`, `mult_sigma1/2/3`, `mult_Ylm_real(l,m,base)`. Overlap solve via `op_DHD.solve<N>`.
- **NOTE:** there is NO `jj_disc_claude.cu` in the tree (the old handoff named one; it does not exist).
  The disc structure must be taken from `jj_corr_block_t`'s disc trace, or written fresh.

## *** RESOLVED 2026-06-12: scheme (A) DIAGONAL (user: "Yes, we use diagonal sum"). ***
Use `meson_pq_wall_v2` for STYLE only; the algorithm is the K-style per-site diagonal estimator below.

### Estimator (free + interacting; valence $D_m=D_{ov}+m$; one VOLUME $Z_2$ source $\eta$ per hit)
Shared: `phi' = D_m^{-1} eta` (one solve/hit; `op_DmH.from_cpu(tmp,eta)` then `op_Dmsq.solve(phi',tmp)`).
`D_m^{-dag} x = (D_m^dag D_m)^{-1} D_m x` -> `op_Dm.from_cpu(tmp,x)` then `op_Dmsq.solve(.,tmp)`.

CONNECTED, channel $a\in\{1,2,3\}$ (s1,s2 spatial, s3 temporal), site $n$, origin $t_0$:
- localize: `src = eta restricted to (n,t0)` (copy eta, zero all but the 2-spinor at site n, slice t0),
  then `src.mult_sigma(a)`  => `src = \sigma_a eta|_{(n,t0)}`.
- `chi_{n,a} = D_m^{-dag} src`  (one solve per (a,n,t0)).
- sink (LOCAL, no solve): `sigphi_a = \sigma_a phi'` (one copy+mult_sigma per a, all sites/times).
- `conn_a(n,t) = chi_{n,a}(t,n)^dag . sigphi_a(t,n)` (2-spinor dot at site n, slice t)
  `=> E = tr[\sigma_a(n,t0) D_m^{-1} \sigma_a(n,t) D_m^{-1}]` (EXACT diagonal-in-$n$; same sign as determ).
- accumulate `Cpp_a[dt] += w_site[n] * conn_a(n,t)`, `dt=(t-t0)%Nt`, `w_site[n]=dual_areas[n]`.
  `write_corr` folds `1/(4pi)`; Vmm=conj (non-parity). MATCHES `corr_deter_local` Vpp directly.

CONNECTED ylm tower ($l=0,1,2$; Eq.4.35/4.36): `Sigma_lm = sum_n A_n Y_lm(n) \sigma_3(n)`, `A_n=dual_areas`:
- source: `srclm = eta restricted to WALL t0`, `.mult_sigma3()`, `.mult_Ylm_real(l,m,base)`
  (mult_Ylm_real multiplies by `dual_areas*Ylm = A_n Y_lm` -- exactly Sigma's weight); `chi_lm = D_m^{-dag} srclm`
  (one solve per (l,m,t0); 9 solves/t0 total over l=0..2).
- sink: `sig3phi_ylm = (sigphi_3) .mult_Ylm_real(l,m,base)`  (= Sigma_lm phi', all sites/times).
- `conn_lm(t) = sum_{n@t} chi_lm(t,n)^dag . sig3phi_ylm(t,n)` (per-slice dot); `Cyl[dt]+=conn_lm`;
  after m-sum `*= 1/(2l+1)`; write_corr folds 1/4pi. Matches jj_local_deter ylm.

DISCONNECTED (separate file; CHEAP, 1 solve/hit -- rides phi' only, no source solves):
- `J_a(t) = sum_n w_site[n] eta(t,n)^dag \sigma_a phi'(t,n) = sum_n w_site[n] eta(t,n)^dag sigphi_a(t,n)` (RAW).
- ylm disc `J_l(t) = sum_n A_n Y_lm(n) eta(t,n)^dag sigma_3 phi'(t,n)` (use sigphi_3 weighted by A_n Ylm).
- write raw `J` (no 1/4pi until the disc-correlator analysis, matching K disc).

### Secondary decisions (DEFAULTS chosen; veto if wrong)
- source = full-VOLUME `fill_z2_source` (matches K; volume needed so phi' reaches all sink times t).
- output = per-hit K-style for resume+jackknife:
  conn `data_<ESNID>/corr_local_stoch_nt0<N>_nhits<H>/corr.<k>.h<h>.h5`,
  disc `data_<ESNID>/disc_local_stoch_nt0<N>_nhits<H>/corr.<k>.h<h>.h5`;
  keys `h0/t0_b/{s1,s2,s3}/{Vpp,Vmm}`, `h0/t0_b/ylm/Vpp/l<l>` (stoch ylm key order), disc `h0/disc/{s1,s2,s3}/J`, `h0/disc/ylm/l<l>/J`.
- VECTOR only (s1,s2,s3 + ylm). NO axial yet (determ loc TODO'd it; add later if wanted).
- free + interacting (`--ens-dir` omit=free); valence mass `--mass-re/--mass-im`.
- npole=11 (match HMC/jj_corr_block_t); no mrhs batching in v1 (loop solves; optimize later if slow).

## (historical) open question that was resolved above
The deterministic loc (Eq.4.29) is **insertion-DIAGONAL**: `sum_n w_n tr[\sigma(n,t0) P \sigma(n,t) P]`
with the **same spatial site $n$** at $t_0$ and $t$ (the current at site $n$ correlated with itself).

`meson_pq_wall_v2` uses a **WALL source**, which estimates the **DOUBLE sum**
`sum_{x@t, y@t0} tr[\sigma S(x,y)^dag \sigma S(x,y)]` (independent sums over source-slice $y$ and
sink-slice $x$) = the $q=0$ projected correlator. The project memory records that **"the $q=0$
double-sum was wrong"** and the chosen observable is the diagonal. So:

- **(A) DIAGONAL (matches `jj_local_deter` / Eq.4.29 exactly):** must use the K-style scheme -- a
  per-site source solve `psi_n(t0)=D^{-dag}\sigma(n,t0)eta`, sink `\sigma(n,t)phi'`, same volume eta both
  ends. Cost ~ `n_sites` source solves/hit (like K tp). Reproduces the deterministic `corr_deter_local_L1`
  and the trio CFT checks ($G_s/G_t\to-1$, ylm tower). This is the physically-validated observable.
- **(B) WALL double-sum (literal meson_pq_wall):** ~1 solve/hit (cheap), but a DIFFERENT ($q=0$)
  observable than the diagonal Eq.4.29; will NOT match `corr_deter_local`.

**Likely intent:** "implement the **local source**" (vs the meson's wall source) = scheme (A), i.e. use
`meson_pq_wall_v2` only as a STYLE template and swap the wall source for the local/diagonal scheme so it
reproduces the deterministic loc. **CONFIRM (A) vs (B) with the user first.** (My recommendation: A.)

Secondary questions to confirm:
1. Source type for scheme A: full-volume `fill_z2_source` (K's choice) vs per-timeslice wall. Volume is
   simplest and matches K; the diagonal-isolation works with either.
2. Output dirs / keys: mirror K's per-hit layout `data_<ESNID>/corr_local_stoch_nt0<N>_nhits<H>/corr.<k>.h<h>.h5`
   with keys `h0/t0_b/{s1,s2,s3}/{Vpp,Vmm}` + `h0/t0_b/ylm/Vpp/l<l>` (stoch ylm key order) + disc
   `h0/disc/...`? Or the determ-style single-file keys? (Propose: per-hit, K-style, for resume + jackknife.)
3. Mass cases: free + interacting `--ens-dir`; valence `--mass-re/--mass-im`. Vector only (no axial yet;
   the determ loc TODO'd axial). Confirm whether axial-loc is wanted now or later.
4. Two files: `jj_local_stoch_claude.cu` (conn, incl. ylm) + `jj_local_disc_stoch_claude.cu` (disc). Or
   conn-file also writes disc (as K does) and the "disc file" is the cheap many-config disc-only path?
   The user said "separate file each" -> propose conn-only file + disc-only file, each standalone.

## Proposed file plan (methodology split; VECTOR currents first)
Three methodology-separated files (table in `project_state.md` 2026-06-12):
- `jj_local_stoch_claude.cu` -- CONNECTED loc **sp/tp only** (s1,s2,s3 diagonal, scheme (A)), per-hit.
  **No ylm here** (ylm is its own file -- user: ylm constrained to local currents, separate file).
- `jj_local_ylm_stoch_claude.cu` -- CONNECTED loc **ylm tower** via the ONE-END TRICK (wall source),
  almost a copy of `meson_pq_wall_v2_claude.cu` (the $Y_{lm}$ double-sum wants the all-pairs structure).
- `jj_local_disc_stoch_claude.cu` -- DISCONNECTED loc (s1,s2,s3 + ylm) loop traces, time+spin dilution,
  almost a copy of `saved_scripts_claude/disc_claude.cu`.
- all `#include "sparse_dirac_claude.h"` (CSR bucketing), `-DN_REFINE_CLI=<L>`.
- run script `tmp_*_claude.sh` (free L=1) + comparison notebook `jj_local_stoch_validate_claude.ipynb`
  (stoch loc vs `corr_deter_local_L1` + CFT), analogous to `jj_free_stoch_validate_claude.ipynb`.

## Implementation chunks (VECTOR, one-by-one)
1. **conn sp/tp file** `jj_local_stoch_claude.cu`: CLI (mirror jj_local_deter + nhits), `w_site=dual_areas`,
   per-hit `eta`; `phi'=D^{-1}eta`; per site $n$ source solve + LOCAL sink `\sigma(n,t)phi'`; build
   `conn(a,t)` for $a=1,2,3$; write per-hit `Vpp/Vmm` + `1/4pi`. **No ylm.** Files: jj_local_stoch_claude.cu.
2. **ylm file** `jj_local_ylm_stoch_claude.cu`: copy `meson_pq_wall_v2`, swap the meson contraction for the
   $g_l$ tower (wall source -> mult_sigma3 -> mult_Ylm_real -> solve -> per-slice contraction). Files: jj_local_ylm_stoch_claude.cu.
3. **disc file** `jj_local_disc_stoch_claude.cu`: copy `disc_claude.cu` (time+spin dilution, all 4 sigma
   channels via accumulate_loop_gamma); adapt weight to `w_site=dual_areas`, add ylm loop `J_l(t)`, jj
   output keys, use $D_m$. Files: jj_local_disc_stoch_claude.cu.
4. **run + validate**: `tmp_claude.sh` free L=1; `jj_local_stoch_validate_claude.ipynb` vs determ + CFT.

## Future TODOs (deferred; marked 2026-06-12)
- **LOCAL AXIAL current** (sp/tp/ylm axial): requires a little care -- the GW factor $(1-D_{ov})$ for the
  axial dressing (jj_local_deter notes it $\to 1$ as $O(a)$, dropped; the stochastic axial-loc needs this
  handled properly). Deferred -- VECTOR currents first, one-by-one. (Exact-K axial already dilutes along in
  `jj_corr_dilute`; this TODO is specifically the LOCAL axial.)
- **LOCAL disc further dilution**: extra dilution (beyond time+spin) can cut the loc-disc noise; possible
  but NOT pursued now.
- (later) interacting ensembles, valence mass (m_F / m_P), reweighting R.

## Status of the prior (K) stochastic validation -- DONE this session
- `jj_free_stoch_validate_claude.ipynb` built + executed: **stochastic K reproduces deterministic
  `corr_deter_exactsum_L1` within errors** (G_t=-1.43e-4+-1.9e-5 vs determ -1.57e-4; G_s=+4.43e-4+-1.6e-5
  vs +4.26e-4; ratio -3.09+-0.41 vs -2.71 @dt5). Signs correct (G_t<0, G_s>0).
- The exact ratio is finite-$t$: sweeps +2.2(dt1) -> -12(dt3 flip) -> -2.7(dt5) -> -1.1(dt12) -> ~-0.9,
  crossing $-1$ near dt~14 (loc is $-1$ at all dt; disp settles ~-1.27 at L=1). dt=5 is NOT asymptotic.
- Data: free stoch L=1 = `corr_nt02_nhits64` (64 hits done); determ `corr_deter_exactsum_L1` COMPLETE;
  exactsum L2 was STOPPED (partial `.tmp`).
- `tmp_claude.sh` currently = the free-stoch L=1 64-hit run (done).
