# Final Production Analysis -- Agent Coordination

Organizer: **qed3-7f** (this session). Reach the organizer as `qed3-7f`.
Started: 2026-09-14.

> TREE LAYOUT (important): the LIVE code + data + blackboard all live in the NESTED inner tree
> `/mnt/barracuda22/qed3/qed3/src/production/` (the outer `/mnt/barracuda22/qed3/src/production/`
> is otherwise empty). This final-analysis dir therefore lives INNER:
> **`/mnt/barracuda22/qed3/qed3/src/production/final/`** -- it holds all scripts, figures, and
> markdowns for the final analysis, sitting next to the code/data it reads.

Blackboard for measurement/ensemble status:
`/mnt/barracuda22/qed3/qed3/src/production/redo_ensembles_claude.txt`.

> DIRECTORY POLICY (NM, 2026-09-14) under `final/`:
> - `peer/` -- ORGANIZER (qed3-7f) files, incl. this coordination doc. Written by 7f only.
> - `shared/` -- WRITABLE by all agents; the common drop-zone (e.g. the shared/common effmass core,
>   cross-channel tables, anything meant for everyone).
> - `analysis_{target}/` -- each analysis agent's OWN dir (scripts, figures, findings, impl plan);
>   written by that agent only.
> - ALL directories are mutually READABLE by every agent. Write only to your own dir + `shared/`.
> This doc now lives at `final/peer/coordination_claude.md`.

> NOTE: the four analysis agents' project memories are ~75 days old and point at the OBSOLETE
> `src/both_3d` tree. The live tree is `src/production`. Each agent must re-verify its scope
> against on-disk data before quoting numbers.

## Roles

### Organizer
- **qed3-7f** -- coordinates the four analysis agents, maintains this doc + the deconfliction map,
  routes cross-channel questions, and relays deliverable direction from NM. Does not itself run the
  per-channel analyses.

### Independent streams (self-managed; coordinate at boundaries only)
- **qed3-6d** -- local MEASUREMENT streams (conn + disc $Y_{lm}$ towers), GPU0 (MPS 2-pack).
  Currently: stride-2 conn fill (classes 3,5,7,9) L3 at0.2 then L1 at0.1
  -> `data_<ens>_vm.../corr_ylm_conn_t00_nhits1_s1/`. Disc remainder (L1 at0.1 + L4 Nf4 g2.0) done.
  HMC generation is done except L4 on SCC. **Hands off GPU0 + its output dirs.**
- **sigma-meson-overlap-sweep** -- two-meson states / distillation, the $\sigma^2$/$0^{++}$
  FOUR-point sector, and the coupled $\{F^2,\sigma^2\}$ mixing GEVP. GPU1 (MPS 2-pack), ~7-day
  L2 perambulator run -> `data_<ens>/distill_Nv24_v2/`. **Hands off GPU1 + that output dir.**

### Analysis agents (report to organizer)
| Agent | Channel | Status |
|-------|---------|--------|
| qed3-d4 (Conn A) | axial two-point (connected only) -- Hankel ell-spectrum | scope in / delivered |
| qed3-60 (rotational symm study) | AXIAL co-worker: m-variational / rotational-symmetry study | new 2026-09-14 |
| qed3-4e | AXIAL co-worker: sp (tangential s1+s2 local-proxy) piece | new 2026-09-14 |
| qed3-a6 (Fin: Fsq) | GLUE co-worker: F^2 / F^4 (0++ scalar glueball) | new 2026-09-14 |
| qed3-42 (jackknifer) | CROSS-CHANNEL aggregator: renormalized dimensions from jk dumps | new 2026-09-15 |
| Fin: Stress tensor | SIGMA-THREAD offshoot (not a final-channel): fermionic T_00 (0++ Delta=3) from distillation, probing the sigma^2 {1,1,1,1} excited 0++; uses the frozen Hankel core, coordinates with sigma-meson-overlap-sweep | new 2026-09-16 |
| qed3-1d (Fin:V)  | vector two-point (conn + disc)   | scope in |
| qed3-a7          | gluonic two-point (glueball / $F^2$ / $F^4$) | scope in |
| qed3-5b (Fin)    | scalar two-point ($\sigma_{PS}$ / $\sigma_{FS}$) | scope in |

## PHASE 2: Renormalized dimensions (space-time anisotropy renorm) -- NM 2026-09-15
GOAL: convert measured dimensionless a_t*m into RENORMALIZED dimensions/masses by renormalizing the
space-time anisotropy. The CONSERVED CURRENTS are the renorm conditions: FERMIONIC sector renormalized
by the AXIAL current (reference = d4 axial tp ell=1 T1); GLUONIC sector by the F field-strength current
(reference = a7 linear F l=1 T1). Fermionic operators (scalar 5b, vector 1d, axial ell>=2 d4, per-m
qed3-60, sp qed3-4e) use the axial current; gluonic operators (F^2/F^4 a6) use the F current. Plus ONE
relative fermion<->gluon dimension, MEASURED DIRECTLY (NM: "we just measure it").
CONSTRUCTION (proposed, jackknifer to CONFIRM exact form + Delta_J values with NM): Delta_O =
(a_t*m_O / a_t*m_{J_s}) * Delta_{J_s}; the a_t cancels in the ratio (= the anisotropy renorm); form the
ratio PER jk sample (aligned by config) then jackknife for the error. Relative dim = m_axial/m_F direct.
AGENT: qed3-42 (jackknifer) aggregates. Correlator agents DUMP per-jk-sample a_t*m arrays + metadata
(L,Nf,gsq,at,binsize,kmin,ell,irrep,op_label,is_current,jk_samples,central,comb) to
final/shared/jk_dumps/<channel>/. Agents still mid-deliverable FINISH first, then dump. GOTCHA: fermion
binsize=10 vs glue binsize=80 -> jk samples NOT aligned across sectors (within-sector ratios fine; cross-
sector relative dim needs common binning or independent-error treatment -- jackknifer's impl-plan call).
First pass: massless, at0.2 first (then at0.1), L1-L4 Nf{2,4,6} where current+operators both FINAL.
Jackknifer dir = final/analysis_renorm/; outputs -> final/shared/renorm_dim/ + peer memo. Dump SCHEMA is
being finalized by the jackknifer -> organizer relays it to the correlator agents.
CONFIRMED binning (all 7 acked ready): FERMIONIC uniform DELETE-1 BLOCK binsize=10 kmin=20 nbin=ncfg//10
(d4 axial=ref, 5b scalar, qed3-60 per-m, 1d vector, qed3-4e sp); GLUE binsize=80 kmin=20 (a7 F=ref, a6 F^2).
SCHEMA REQ (5b+a6): dump must carry a per-bin config-id map so operator/current ratios are provably on the
SAME jk decomposition (same configs + bin assignment); jackknifer verifies alignment + ensemble coverage
before any ratio. Ready-now: d4/5b/a7/a6/qed3-60 (at0.2). Finish-first: 1d (repoint vector driver), qed3-4e
(sp on NM hold). d4 re-dumps after 6d conn topup.
SCHEMA jk_dump_v1 FINAL 2026-09-15 (jackknifer): spec `final/analysis_renorm/jk_dump_schema_claude.md`
+ validator `jk_dump_validate_claude.py` (run before pinging) + impl plan `renorm_dim_impl_plan_claude.md`.
File per (channel,at): `final/shared/jk_dumps/<channel>/<channel>_at<at>_claude.json`. Channels: axial_tp(d4),
axial_perm(60), axial_sp(4e), scalar_ps(5b), vector(1d), glue_F(a7 l=1), glue_Fl2(a7), glue_Fsq(a6).
CONTENT = per-jk FITTED plateau a_t*m (delete-1, published estimator/window/weights) + jk_samples_sys
(+1-window shift) + central_sys + MANDATORY cfg_k (all used ids asc) + bin_k (per-bin config-id list,
len==njk). Ratio accepted ONLY when bin_k exactly equal, else independent-error fallback (flagged).
Fermion = odd-k stride-2, count-bin10 on sorted(glob) LEXICOGRAPHIC (d4/5b/60 align); glue all-k count-bin80
numeric (a7/a6 align, Nc=3980/nbins=49). CROSS-sector rho = m_axial(ell1)/m_F(l1) starts INDEPENDENT-error
(count-bins never align); possible common k-interval binning (W=80) addendum for d4+a7 pending NM.
OPEN w/ NM (jackknifer, non-blocking): pure-ratio form, Delta_J=2 both currents, exact rho estimator, binning policy.
ADDENDUM v1.1 (NM 2026-09-15, RESOLVED): UNIFIED bin size 80 as a k-INTERVAL throughout. bin b = configs with
20+80b <= k < 100+80b (NUMERIC k, kmin=20); drop last bin only if incomplete; sample b = delete all of bin b;
bin_rule={type:kinterval,W:80,order:numeric}. FERMION agents RE-DUMP with it (d4 axial_tp both at [still the
critical-path reference], 5b scalar_ps both at, qed3-60 axial_perm, 1d vector [numeric now, not lexicographic],
4e axial_sp when settled); fermion conn = odd-k only => ~40 cfg/bin; re-dump IN PLACE + bump dump_date;
dump central/stat shift slightly vs native bin10 (published mds unchanged). GLUE = NO re-dump (jackknifer
verified all 114 glue entries already sit on W=80 k-slabs). DEFINITION (NM-set): pure ratio
Delta_O = (m_O/m_J)*2, Delta_J = 2 for BOTH currents; rho = m_axial(ell1)/m_F(l1) per ensemble, correlated jk
once fermion re-dumps land. Scope = whatever each agent delivers FINAL. Schema/validator updated in
final/analysis_renorm/.

## Per-channel scope (summary)

### qed3-d4 -- Axial two-point (CONNECTED only)
- Connected axial current only. The disc driver computes the VECTOR current only, so the axial
  sector is purely connected -- no overlap with qed3-6d's disc work.
- tp = s3 ($J\cdot\hat n$): faithful scalar, ell=1(T1), 2(H), 3(T2+G). sp = s1+s2: local proxy, ell=0,1 (both $\Delta=2$).
- Ensembles: L1-L4 sphere, Nf 2/4/6, gsq scan, massless + Family-B massive
  (mRe $\in\{0.0106,0.0529,0.1057,0.2114\}$); at0.2 & at0.1 paired (L1 g1.0, L2 g2.0). Cut ncfg>=100.
- Inputs: `corr_ylm_conn_t00_nhits1_s1` h5 + `axial_tp_masses_summary_claude.md` (single source of truth).

### qed3-1d -- Vector two-point (conn + disc)
- Local-current $Y_{lm}$ tower. tp = s3 ($J\cdot\hat n$): true scalar, faithful, reliable.
  sp = s1s1+s2s2: local proxy only (VSH mixing; needs covariant link current $G^s$ for a real sp).
- Physical correlator = $-C_\text{conn} + C_\text{disc}$; conn-only = $-C_\text{conn}$.
  Disc uses tb=2 dilution, DC-subtracted (LSD Eq.12) + large-t plateau ($16\le dt\le Nt/2$) subtracted.
- Must confirm which vector conn+disc data exist in the production tree (redo campaign is L1-L4).

### qed3-a7 -- Gluonic two-point (glueball / $F^2$ / $F^4$)
- Icosahedral-orbit spatial Wilson-loop SHAPES basis (`includes/wilson_shapes_claude.h`), shapes 0-6,
  face_sign ON. Channels: linear $F_{12}$ (pseudoscalar l=1 T1) + $F^2$ ($0^{++}$ l=0);
  $F^4$ added (op_pow p={2,4} BEFORE the $Y_{lm}$ sum). Driver `glue_f2_v2_shapes_claude.cu`.
- GEVP via `glue_gevp_analysis_claude.o` (host streaming-jk). Prod nops2=4 vacsub=0 => $0^{++}$ is STATE 2.
- Free-limit anchors (validated): $F$ l=1 = $\sqrt2$; $F^2$ $0^{++}$ = $2\sqrt2$ (two-photon).

### qed3-5b -- Scalar two-point ($\sigma_{PS}$ / $\sigma_{FS}$)
- $\sigma_{PS}$ conn $V_{++}=\text{tr}[W_0 D_m^{-1} W_t D_m^{-1}]$; $\sigma_{FS}$ furnished (both $(1-D_\text{ov})$
  insertions) $V_{--}^{FS}=\text{tr}[W_0 \Gamma W_t \Gamma]$, $\Gamma=(1-D_\text{ov}^\dagger)D_m^{-\dagger}$.
  Scalar-density $Y_{lm}$ tower (a=0 identity vertex). Disc DROPPED ($\langle\sigma\rangle=0$ at m=0).
- Drivers: `jj_local_ylm_scalar_conn_stoch_claude.cu`; analysis `jj_scalar_ylm_analysis_claude.ipynb`.

## Deconfliction map
- GLUE co-worker split (a7 + qed3-a6, share analysis_gluonic/, handoff_fsq_claude.md DONE): a7 keeps
  linear F_12 (l=1) and stays knowledgeable on it; qed3-a6 owns F^2 / F^4 (0++, l=0, glue_f2_v2_shapes).
  Boundary CONFIRMED by NM 2026-09-14: a7 ALSO keeps Fl2 (l=2) -- same linear-F operator glue_msm_shapes
  (lsel 1 vs 2). So a7 = linear F l=1 + Fl2 l=2; a6 = F^2 / F^4 (0++). Distinct output prefixes (a7
  gevp_Ffix_*/F_fixt0_*; a6 *Fsq*); shared binary/driver signatures kept stable.
- AXIAL co-worker split (d4 + qed3-60 + qed3-4e, share analysis_axial/): qed3-60 owns per-m /
  rotational-symmetry (effmass_axial_tp_l3_perm*, per-m block-Hankel, per-m tuning; handoff_mvariational_claude.md);
  qed3-4e owns the sp (tangential s1+s2 local-proxy) piece (handoff_sp_claude.md, mid-study); d4 keeps the
  m-AVERAGED Hankel ell-spectrum TP + fit/summary (offloaded sp so she stays tp-focused w/o compactifying).
  HANDOFF DONE 2026-09-14: qed3-4e has her own forked sp driver (hankel_ell_spectrum_sp_claude.py,
  default --channel sp, own cache -> fully decoupled from d4's tp driver). Each of the 3 has own driver +
  own scratchpad cache. ONLY shared file = Hankel CORE effmass_axial_tp_l3_perm_hankel_claude.py (H),
  now FROZEN -- all three coordinate signature-level edits.
  sp ELL RULE (d4 clarified): sp KEEPS ell=0,1,2 -- the tp ell=0 conserved-charge drop does NOT apply to
  sp (sp ell=0 is a real tangential Delta=2 state). sp correlator is sign-flipped; VSH-mix near-degeneracy
  expected; sp params start from tp per-L ladders/windows verbatim (not yet sp-tuned).
- **Coupled $\{F^2,\sigma^2\}$ mixing GEVP -> owned by sigma-meson-overlap-sweep.** Glue (a7) and
  scalar (5b) supply single-channel two-point inputs only; they do NOT run the coupled mixing GEVP.
  $F^2$ mixes only with parity-EVEN $\sigma^2$ (P+ {PP,FF}); $\langle F^2\,FP\rangle=0$ by parity.
- **a_t artifact + valence a_t bug fix** (`project_valence_at_fix`): shared caveat across axial (d4),
  vector (1d), scalar (5b) fermion masses at at0.1. The at0.1 fermion points used a WRONG valence
  a_t=0.2 (driver hardcode); recompute in progress. Re-verify any a_t-artifact conclusion after the
  corrected at0.1 data lands. Coordinate the re-verify across the three fermion channels.
  CONFIRMED 2026-09-14 (qed3-5b): the at=0.2 hardcode is present in the SCALAR analysis scripts too
  (effmass_scalar_expfit_claude.py, scalar_over_axial_vs_gsq_claude.py build only at0.200000 names)
  => current scalar mass tables are at0.2-ONLY, even though at0.1 scalar data exist on disk (36 at0.1
  dirs, scalar keys confirmed). Same class of bug in the fitters, not just the measurement driver.
- **tp = s3 ($J\cdot\hat n$)** is the faithful scalar proxy for BOTH vector (1d) and axial (d4);
  sp is local-proxy-only in both. Keep this convention consistent across the two channels.
- **Scalar 2pt (5b) vs $\sigma$-meson 4pt/distillation (sigma agent):** 5b stays strictly in the
  two-point sector; the four-point / $\sigma^2$ / $0^{++}$ mixing is the sigma agent's.

## Canonical a_t-hardcode fix pattern (from qed3-6d, territory owner)
Two hardcodes per fitter; the SECOND is the subtle correctness one:
1. ENS-DIR glob: literal `at0.200000` -> `at%.6f` with `at` a loop variable (else at0.1 dirs never
   built => at0.1 data silently UNREACHABLE / skipped).
2. The `/at` physical-mass conversion: `me = meff_acosh(C) / at` uses a MODULE-LEVEL `at=0.2`.
   Must be the PER-ENSEMBLE a_t, else at0.1 physical masses are off by 0.2/0.1 = 2x. This is the
   analysis-side analog of the driver kappa_t bug. (physical m = acosh_effmass / a_t.)

BROKEN exemplars:
- scalar (5b): effmass_scalar_expfit_claude.py (:18 at=0.2, :39 at0.200000, :79 /at);
  scalar_over_axial_vs_gsq_claude.py (:15,:50,:90).
- vector (1d): effmass_vec_expfit_claude.py (:16 at, :34-35 esn literal, :111 /at, :31 tt=dtp*at);
  jj_ylm_prod_analysis_redo_claude.ipynb (hand-edited nb -> surgical/coordinate).
- axial (d4, flagged by 1d): effmass_conn_claude.py + effmass_plateau_tables_claude.py (same shape).

REFERENCE scripts already CORRECT -- copy verbatim: effmass_axial_tp_at_claude.py
(:23 `ATS=[0.2,0.1]`, :28 `at%.6f`, :64 `/at` with loop at), ratio_vs_at2_claude.py,
mat_vs_scales_claude.py, effmass_axial_tp_at_{expfit,ratio}_claude.py, effmass_l2nf6_at_claude.py.
Pattern: `for at in ATS`, build dir with `at%.6f`, divide by that loop `at`.
Non-loop parse: `at = float(re.search(r'at(0\.\d+)', ensdir).group(1))`.

at0.1 EXISTS ONLY for: L1 Nf{2,4,6} x g{0.5,1.0,1.5} and L2 Nf{2,4,6} x g2.0. NO L3/L4 at0.1.
PER-a_t FIT WINDOWS DIFFER (don't reuse one window): axial at0.1 dt[1,2.4] vs at0.2 dt[1.4,4.0].
Scalar/vector fixers must set at0.1 windows too, not inherit at0.2.

DATA STATUS (corrected-a_t): at0.1 conn k=1 mod10 complete everywhere, classes 3,5,7,9 filling now on
GPU0; at0.1 disc L1 (both hits) + L2 done. OLD wrong-a_t at0.1 fermion h5 ARCHIVED to
<data_dir>/archive_at0.2/; live at0.1 dirs hold corrected-a_t data (driver fix commit bca707b:
at_from_ensdir + --at in jj_local_ylm_scalar_conn_stoch{,_fnal} + jj_local_ylm_disc_stoch). Condensate
driver still UNFIXED = out of scope. Refs: project_valence_at_fix + redo_ensembles_claude.txt.

## a_t-hardcode audit results (2026-09-14)
KEY REFRAMING: in the fermion analysis fitters the at=0.2 hardcode is a clean OMISSION, not a
mis-scaling -- each broken fitter globs the literal at0.200000 dir AND divides by at=0.2, so name-a_t
== divisor-a_t (SELF-CONSISTENT); at0.1 is simply SKIPPED (never loaded), never mislabeled. So the
existing at0.2 tables are correct for their at0.2 scope; the gap is MISSING at0.1 coverage.
=> Adding at0.1 coverage for scalar/vector = NEW work (create at-parameterized fitters from the axial
_at_ reference pattern), which is a DELIVERABLE-SCOPE decision for NM, not an automatic fix.

- AXIAL (d4): fitters need NO edit. 8 at0.2-only production fitters skip at0.1 cleanly; at0.1 is
  already handled by 3 CORRECT at-parameterized scripts (effmass_axial_tp_at{,_expfit,_ratio},
  ATS=[0.2,0.1]). No kappa_t logic in any .py. REAL exposure: the at0.1 axial DATA those 3 scripts
  read may carry the driver valence-a_t=0.2 bug -> axial a_t-artifact conclusion suspect until 6d
  recomputes. Axial action = downstream RE-RUN after recompute, no code fix.
- VECTOR (1d): 2 broken fitters -- effmass_vec_expfit_claude.py, jj_ylm_prod_analysis_redo nb. NO
  at-parameterized vector fitter exists => at0.1 vector coverage needs a NEW _at_ fitter (copy axial).
- SCALAR (5b): 5 scripts -- effmass_scalar_expfit (primary), effmass_scalar_L3, scalar_over_axial_vs_gsq,
  + cross-channel effmass_conn & effmass_full. NO at-parameterized scalar fitter => same as vector.
  12 at0.1 scalar ensembles on disk (L1 Nf{2,4,6}xg{0.5,1,1.5} + L2 Nf{2,4,6}xg2.0, 474-1000 cfg).
- GLUE (a7): audit DONE. Mass conversion CLEAN -- glue_gevp_analysis takes `at` as CLI arg; at0.1
  handled by run_glue_gevp_at01_claude.sh (AT=0.1). GAP: at0.1 conversion script exists only for LINEAR
  F (glue_msm_shapes), NOT for F^2/F^4 (glue_f2_v2_shapes) -> need an AT=0.1 F^2 pass if at0.1 F^2 wanted.
  Ensemble-name hardcode bypassed in prod (explicit ens_dir arg7). GLUE-SPECIFIC driver issue (6d
  territory): glue_f2_v2_shapes_claude.cu:254 `const double at=0.2` feeds beta_s = at/(vol*gsq)
  (LINEAR in at, action_ext:366) which drives the spatial Wilson FLOW => at0.1 glue h5 were flowed with
  beta_s(0.2) = 2x the at0.1-consistent coupling (same absolute smearing for all a_t). Latent smearing
  INCONSISTENCY at the smearing level.
  CORRECTED 2026-09-14 (a7 + NM, RESOLVED): this does NOT caveat the at0.1 glueball MASSES. The driver
  `at` enters ONLY the spatial per-timeslice Wilson-flow smearing (beta_s propto at); smearing changes
  interpolator overlap / signal-to-noise, NOT the temporal transfer matrix, so GEVP/plateau MASSES are
  smearing-level-INDEPENDENT -- existing at=0.2-flowed at0.1 h5 give the SAME masses as an at=0.1-flowed
  re-dump. The only genuinely at-dependent result step is the mass conversion Delta_eff=-log(lam)/(dt*at),
  which uses the ANALYSIS at (CLI AT=0.1) and was already correct. So existing at0.1 glue is NOT caveated
  for masses; the _v2.1 re-dump was shelved as unnecessary (NM's call). a7 offers an optional cheap 30-cfg
  empirical confirm (at=0.2 vs at=0.1 flow -> plateaus agree); NOT a blocker. a7 still fixed the glue
  flow-at across all 3 shape drivers (at_from_ensdir; at0.2 bit-identical) -- available, harmless.
  Name-builder hardcode harmless in all three (explicit ens_dir arg bypass).

ALL FOUR AUDITS COMPLETE 2026-09-14. Fermion side: no active mis-scaling; at0.1 is a clean omission
(new _at_ fitters needed for scalar/vector IF at0.1 in scope). Glue side: mass conv clean, one latent
flow-smearing inconsistency at at0.1 + F^2/F^4 at0.1 conversion-script gap.

OWNERSHIP (settled with 6d, all GATED on NM confirming at0.1 in scope):
- Fermion at-enable: 6d owns the SHARED cross-channel fitters effmass_conn + effmass_full (single edit,
  avoids d4/1d/5b triple-editing). Per-channel non-shared _at_ fitters: the owning analysis agent.
- Glue at0.1 flow fix + diagnostic: a7 owns (glue shape drivers are the glue thread's, NOT 6d's
  fermion-only project_valence_at_fix). Spans all 3 drivers (v2 :254, glue_f2_shapes :231, glue2_msm :241).

DATA STATUS (6d, 2026-09-14): at0.1 fermion CONN usable NOW at good stats (~600-800/1000 cfg per L1
at0.1 ens: k=1 mod10 =200 complete + classes 3,5,7,9 =400-600 partial from FNAL rerun); at0.1 DISC done
(L1 both hits, L2). Full stride-2 topup to ~1000 = PHASE B, runs AFTER PHASE A (L3 at0.2); L3 is long
pole => ~2-3 weeks out. Runner phase order (L3-first vs at0.1-first) is flippable = NM's call.
- CROSS-CHANNEL scripts (shared, need coordinated ownership -- flagged to 6d):
  effmass_conn_claude.py (axial tp l=1,2 + scalar PS/FS l=0,1) and effmass_full_claude.py
  (vector/axial/scalar). A fix here touches all three fermion channels at once.
- NOTEBOOKS (6d territory / hand-edited, surgical only): jj_ylm_prod_analysis_redo_claude.ipynb,
  conn_massless_redo_analysis_claude.ipynb.

## Scope decisions from NM (living)
- 2026-09-14: FINAL analysis is MAINLY MASSLESS. The Family-B massive (mRe) axis is DEFERRED; NM will
  direct when the exceptional massive ensembles are done. => all four channels focus massless for now;
  drop massive-mass channels unless NM says otherwise. (Vector already massless-only on disk.)
- 2026-09-14: at0.1 scoping is LEFT OPEN, to clarify as the work proceeds. Default working: at0.2 first.
- 2026-09-14 (UPDATE): at0.1 is now IN SCOPE (NM lifted the hold, "go ahead"). Valence-propagator
  correctness settled (6d recompute, commit bca707b; wrong-a_t archived). => the at0.1 a_t-artifact arm
  is greenlit: fermion at0.1 grid (L1 g{0.5,1,1.5} + L2 g2.0, Nf2/4/6, all 12 OK 468-990 cfg).
  Fermion fitter at-enable (6d owns shared effmass_conn/effmass_full) is un-gated; per-channel _at_
  fitters proceed under the common method. GLUE at0.1: the beta_s(0.2) flow-smearing does NOT gate the
  glueball MASSES (CORRECTED 2026-09-14, a7+NM -- see the glue-audit correction note above); the ONE
  real glue at0.1 item is the 4-ensemble MEASUREMENT gap (a7 prepared a user-run handoff).
- 2026-09-14 (6d, on-disk verified): at0.1 CONN DATA-READY confirmed. All 12 at0.1 conn ensembles'
  live corr_ylm_conn_t00_nhits1_s1 at at0.100000 hold ONLY corrected-a_t h5 (oldest mtime post commit
  bca707b 2026-08-21 14:31; wrong-a_t archived to archive_at0.2/ for the 6 that had it; g0.5/g1.5 L1
  first-measured post-fix so never wrong). disc at0.1 corrected + done (L1 both hits, L2). STATS CAVEAT
  (not correctness): k=1 mod10 complete (200/ens); classes 3,5,7,9 PARTIAL for L1 (400-600/ens, topping
  up PHASE B after L3), ~complete for L2 (~800/ens) => values trustworthy now, at0.1 error bars shrink
  as 3,5,7,9 fills. BOTH greenlights (NM scope + 6d data) IN => d4 released for at0.1 axial.
- 6d owns the shared cross-channel fitter at-enable (effmass_conn + effmass_full) and will confirm the
  edit-go directly with NM before editing (does not act on a relayed greenlight as edit-authorization).
- 2026-09-14 (NM directive, relayed via d4): OMIT ell=0 from the VECTOR channel spectrum. Physics: for
  a conserved vector current the ell=0 tp component is the total charge -- t-independent/non-propagating,
  no mass. => 1d drops ell=0 from the vector tower. RESOLVED 2026-09-14: NM confirms AXIAL-vector ALSO
  excludes ell=0 (same conserved-current charge/non-propagating reason). => ell=0 OMITTED for BOTH the
  vector and axial CURRENT towers (axial now shows ell=1,2,3; at0.2 re-emitted 108 rows, at0.1 re-emit
  after its run). SCOPE NOTE: this omission is CURRENT-TOWER-SPECIFIC. Scalar (sigma_PS/sigma_FS) and
  glue (F^2 0++) are NOT conserved currents -- their ell=0 IS the physical ground state and MUST be kept.
  Do NOT drop ell=0 from scalar or glue. REFINEMENT (d4, 2026-09-14): within AXIAL, the ell=0 drop is
  TP-SPECIFIC (the J.nhat conserved-charge monopole) -- the SP (tangential) piece KEEPS ell=0,1,2 (sp
  ell=0 is a real tangential Delta=2 state, not the conserved charge).

## Conventions
- Each agent's `*_impl_plan_claude.md` (and its scripts, figures, findings docs) lives INSIDE that
  agent's own analysis directory `final/analysis_{target}/` -- not at the final/ root, not elsewhere.
- THERMALIZATION CUT (kmin) = 20 THROUGHOUT (NM, 2026-09-14; CORRECTED from an initial 100): discard
  configs with k < 20 (first-20 thermalization) for every stream, all channels/ensembles. This is
  SEPARATE from the ncfg>=100 minimum post-cut config count -- do not conflate the two.
- OUTPUT IN LATTICE DIMENSIONLESS UNITS (NM, 2026-09-14, via d4): report masses/spectra as the
  dimensionless lattice value (e.g. Hankel/GEVP: log(lam_t/lam_{t+1}), i.e. a_t*m, NO /a_t division);
  convert to physical later by /a_t. Rationale: the a_t-artifact then reads directly -- the dimensionless
  a_t*m is ~a_t-INDEPENDENT (flat) if it is an artifact, and DOUBLES in physical at half a_t. So at0.1
  vs at0.2 are compared in dimensionless a_t*m. Applies to all channels for cross-a_t comparison.
  UNITS-LABEL TRAP (5b, likely all fermion channels): the acosh effective mass is ALREADY a_t*m, so old
  fitters that do `meff_acosh(...) / at` (e.g. effmass_scalar_expfit:79) actually emit PHYSICAL m while
  their header says "a_t m0" -- a MISLABEL under the new convention. Under the common method report the
  RAW acosh (a_t*m), drop the /at, and DO NOT trust old "a_t m0"-labeled tables. Reinforces: regenerate,
  don't trust pre-existing tables.
  PER-CHANNEL DIMENSIONLESS RECIPE: fermion (d4/1d/5b) = raw acosh effmass, no /at. Glue (a7) = run the
  GEVP with AT=1.0 (a_t*m = a_t*Delta_eff = -log(lam)/dt, dt=1); existing gevp_f2_*.dat used AT=0.2 =
  PHYSICAL, so dimensionless = x0.2 -- a7 re-emits baseline with AT=1.0 for the final product.
  RULE (d4): divide by a_t EXACTLY ONCE, at the very END, to get physical. The shared reference method
  effmass_generic_claude.py is FIXED (effmass_jk returns raw acosh = a_t*m, no /at; label = LATTICE a_t*m)
  -- so the fix propagates with the common method. CONFIRMED mislabel in old tables of ALL 3 fermion
  channels (scalar effmass_scalar_expfit:79, vector effmass_vec_expfit:111, axial effmass_generic +
  axial_tp_masses_{generic,summary}): those "a_t m0"-labeled tables are actually PHYSICAL m -- regenerate
  dimensionless, do not trust the labels.

## Completeness tables MUST include at0.1 rows (NM, 2026-09-14)
Each channel's completeness table includes BOTH at0.2 AND at0.1 ensembles (run the shared scanner with
--at 0.1 for the at0.1 grid, merge rows; the `at` column distinguishes them). This is INVENTORY only --
it does NOT greenlight the at0.1 analysis arm (that scope decision is still open); it just makes the
at0.1 config counts visible in the list. at0.1 grid: FERMION = L1 Nf{2,4,6} g{0.5,1.0,1.5} + L2 Nf{2,4,6}
g2.0; GLUE = L1 Nf{2,4,6} g{0.5,1,1.5} (a7 to confirm any L2 at0.1 glue on disk). Vector disc at0.1
exists for L1 + L2 (6d: L1 both hits + L2 done).

### at0.1 rows merged -- ALL 4 tables updated 2026-09-14 (INVENTORY ONLY)
- FERMION at0.1 = 12/12 OK (468-990 cfg after kmin=20): axial, vector (conn+disc), scalar all clear
  ncfg>=100 for the full at0.1 grid (L1 Nf{2,4,6} g{0.5,1,1.5} + L2 Nf{2,4,6} g2.0). Combined tables
  = 48 ens (36 at0.2 + 12 at0.1); the only LOW are the 2 at0.2 L3 points (Nf4 g1.5, Nf6 g4.5 = 98).
  Vector at0.1 notably clears conn+disc everywhere (disc 198-420), unlike at0.2 L3/L4.
- GLUE at0.1 = PATCHIER (a7): only 8 OK (L1 {Nf2 g0.5/1/1.5, Nf4 g1, Nf6 g1} + L2 {Nf2/4/6 g2});
  4 EMPTY = stale <20-cfg partials that fail the kmin=20 cut (L1 Nf4 g0.5, Nf6 g0.5, Nf4 g1.5, Nf6 g1.5);
  rest NO_ENSDIR. GAP: those 4 L1 glue at0.1 ensembles have gauge configs (fermion measured them fine =
  OK) but the GLUE MEASUREMENT is missing => if the at0.1 arm is greenlit, glue needs measurement runs
  on those 4 (a measurement gap; NOT a gauge gap). a7 prepared a user-run handoff for them
  (tmp_glue_at01_gap_claude.sh, both dumpers, at=0.1 flow, resume-safe), flagged to NM as measurement,
  to run after the L4 at0.2 topup (avoid CPU oversubscription). [The beta_s(0.2) flow-smearing does NOT
  caveat glue MASSES -- see the CORRECTED note in the audit section.]

## Completeness (checkpoint-count) table -- COMMON style (d4, propagate verbatim)
Reusable scanner (channel-agnostic, glob-only, no h5 reads): `final/shared/completeness_scan_claude.py`
(authored by d4 in analysis_axial/, copied to shared/ by 7f). Each channel runs it with its own subdirs
and saves the resulting table in ITS OWN analysis_{target}/ dir.
- Invocation: `--conn-subdir <dir>` (+ `--disc-subdir <dir>` for conn+disc channels), `--label` to taste.
  axial/scalar = conn-only; vector = conn+disc; glue = point --conn-subdir at its density/shape subdir.
- Columns: `L | Nf | gsq | at | ncfg_conn | ncfg_disc | status` (conn-only prints n/a in ncfg_disc).
- ncfg = per-CONFIG count (NOT h5 file count -- watch multiplicity, e.g. vector disc=2 files/config)
  AFTER kmin=20 thermalization cut. SEPARATE from the ncfg>=100 minimum that drives the LOW flag.
- status on the RELEVANT count -- conn-only: ncfg_conn; conn+disc: min(conn,disc) (matched analysis is
  conn-disc limited): OK(>=100) / LOW(<100) / EMPTY / NO_CONNDIR / NO_ENSDIR.
- Ordering: L asc, then gsq asc, then Nf asc. Two header comment lines (channel+at+kmin; status legend)
  + a roll-up line at the end.
- Reference table: final/analysis_axial/axial_conn_completeness_claude.md.

COMPLETENESS TABLES -- ALL 4 CHANNELS DONE (at0.2 massless grid, kmin=20), 2026-09-14:
- axial (conn): analysis_axial/axial_conn_completeness_claude.md -- 36 ens, 2 LOW (L3 Nf4 g1.5=98, Nf6 g4.5=98).
- scalar (conn): analysis_scalar/scalar_conn_completeness_claude.md -- 36 ens, same 2 LOW.
- vector (conn+disc): analysis_vector/vector_conndisc_completeness_claude.md -- 36 ens, 18 LOW: L1+L2 OK,
  L3+L4 all disc-limited LOW (disc 35-98). Conn-only clears 100 at L3/L4 except the same 2 marginal 98.
- glue: analysis_gluonic/glue_completeness_claude.md -- 36 ens, ALL OK. L4 at0.2 TOPUP COMPLETE
  2026-09-14 (both dumpers glue_f2_v2_shapes + glue_msm_shapes = 799, Nf6 g6=778; post-kmin=20: 780,
  Nf6 g6=759). Every at0.2 massless glue ensemble L1-L4 now OK for BOTH observables. (Nf6 g6 stream
  still generating -> a7 will resume tmp_claude.sh later to catch newer configs.)
  PENDING (NM-run): the 4-ensemble at0.1 L1 glue MEASUREMENT gap handoff tmp_glue_at01_gap_claude.sh --
  now UNBLOCKED (it was waiting behind the L4 topup, which is done). Ready for NM to launch.
- SCANNER: final/shared/completeness_scan_claude.py now covers ALL 4 channels (a7 patch, additive knobs;
  fermion defaults byte-for-byte unchanged). New flags: --bare (bare gauge dir, no _vm suffix),
  --conn-subdir (now defaults ""), --file-prefix/--file-suffix (default corr / .h0.h5). Glue invocation:
  `--bare --conn-subdir "" --file-prefix glue_f2_v2_shapes --file-suffix .h5 --label glue`. a7 verified
  glue via patched tool = identical to its standalone. VERIFICATION COMPLETE 2026-09-14: the SHARED copy's
  code paths are byte-identical-verified -- conn+disc via 1d (vector), conn-only via 5b (scalar),
  glue-bare via a7. (Axial was regenerated via d4's OWN reference copy, schema-identical; its conn-only
  path == scalar's, already covered by 5b's shared-copy run, so no separate axial-through-shared test
  needed.) => final/shared/completeness_scan_claude.py is the SINGLE completeness tool for all 4
  channels; a7's standalone completeness_scan_glue_claude.py can be retired.
- CONFORMANCE (7f check, NM-requested): vector/scalar/glue tables MATCH the common style. AXIAL saved
  table was OLD schema (L|gsq|Nf|n_tot|ncfg(kmin>=20)|status, no `at`, no conn/disc split) -> d4 to
  REGENERATE via the patched shared tool (fixes schema + doubles as fermion-side patch verification).

## Workflow / method-sharing plan (NM, 2026-09-14)
- START POINT: AXIAL CONNECTED, at0.2, L1-L4. Some conn data are incomplete -- that is expected; build a
  GENERIC effective-mass analysis script that runs on partial stats, not a per-ensemble one-off.
- The effmass METHOD must be SHARABLE (generic across channels). NM does focused work with the dedicated
  axial agent (qed3-d4) to develop it; the resulting method is then shared with the organizer (qed3-7f),
  who propagates the reusable core to the other channels (vector 1d, scalar 5b, glue a7 where applicable).
- Consequence: other analysis agents HOLD on building their own effmass fitters -- adopt the shared
  method when it is propagated, to keep one common approach rather than divergent per-channel scripts.
- PROPAGATION MECHANICS (d4): shared effmass core is authored in analysis_axial/; the organizer (7f)
  propagates by COPYING the reusable core into each channel's own analysis dir (per-channel copy +
  thin per-channel driver). NO cross-dir imports -- each analysis_{target}/ stays self-contained.
- GENERAL PRINCIPLE (5b heads-up, likely all channels): current on-disk ncfg now EXCEED the older
  mass-table counts (data grew since those tables were generated) => REGENERATE all tables under the
  common method; do NOT trust the pre-existing mass tables as-is.
- CHECKPOINT-COUNTING GOTCHA (1d, 2026-09-14): count DISTINCT configs (k>=20), NOT h5 files -- some
  streams write >1 file/config. Vector DISC = 2 files/config (hit index h0,h1) => file count is ~2x the
  config count. Conn (nhits1) = 1 file/config. Apply per-stream multiplicity when building ncfg tables.
- VECTOR ncfg-cut consequence (1d corrected): PHYSICAL (conn+disc) vector correlator meets ncfg>=100
  only at L1+L2 (conn 990-1990, disc 198-1194). At L3/L4 disc ncfg=35-98 (<100) => only CONN-ONLY
  (-C_conn) usable there (conn mostly >=100; marginal 98: Nf4 L3 g1.5, Nf6 L3 g4.5).

## Standing guidance to analysis agents (pending NM's deliverable definition)
- OK now: read-only on-disk inventory to re-verify scope vs `src/production`.
- HOLD until NM sets the deliverable: committing final fit windows; launching NEW measurement runs
  (new runs are handoff scripts to NM, not agent-run).

## Open questions for NM
1. Deliverable/structure of the final analysis (paper tables/figures? unified spectrum? per-channel
   summary docs feeding a combined document?).
2. qed3-5b: cover the condensate (one-point) order params $\sigma_{PS}/\sigma_{FS}$ too, or strictly two-point?
3. In-scope ensembles for the FINAL product (all L1-L4, or a subset) and confidence/priority ordering
   across the four channels.

## Log
- 2026-09-14: organizer qed3-7f appointed. Greeted qed3-6d, sigma-meson-overlap-sweep, and the four
  analysis agents. Collected all four scopes. Created `src/production/final/`. Awaiting NM's next command.
