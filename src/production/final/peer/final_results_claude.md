# Final Production Analysis -- Recorded Final Numbers (organizer memo)

Maintained by organizer qed3-7f. As each channel reports its FINAL dataset, the final-numbers
markdown is recorded/pointed-to here and the final PNGs are copied into `final/shared/`.
Convention: lattice dimensionless a_t*m (divide by a_t once at the very end for physical); kmin=20;
ncfg>=100 LOW flag; ell=0 dropped for vector/axial current towers, kept for scalar/glue.

## RENORMALIZED DIMENSIONS (Phase 2, jackknifer qed3-42)
Definition (NM): Delta_O = 2 * m_O / m_J (pure ratio; Delta_J = 2 for BOTH currents). m_J = axial ell=1 T1
(fermionic sector) / F l=1 T1 (gluonic sector). rho = m_axial(ell1)/m_F(l1) ties the two sectors. Config-
aligned correlated jackknife on the unified W=80 k-interval binning.
FULL at0.2 SET DONE 2026-09-15 (jackknifer; tables in analysis_renorm/, copied to final/shared/renorm_dim/).
Tables: renorm_dim_masses_at0.2_claude.md (324 operator rows -- axial ell2/3 72, scalar PS ell0-3 144,
vector ell1/2 36, per-m ell3 18, glue Fsq 18 + Fl2 36; PLUS rho table 36 ens) + renorm_dim_masses_at0.1_claude.md
(2026-09-16: 72 operator rows -- axial ell2/3 24, scalar PS ell0-3 36, glue Fsq 12; PLUS rho on all 12
at0.1 ens; vector at0.1 still pending 1d. ALL rows config-aligned).
**KEY RESULT (jackknifer, at0.1 vs at0.2): the renormalized Delta is a_t-STABLE within a few % -- i.e. the
conserved-current anisotropy renormalization REMOVES the bare a_t-artifact.** At L1: scalar ell0 1.84-1.94
(at0.1) vs 1.87-1.96 (at0.2); scalar ell1 ~2.53 vs 2.47-2.62; axial ell2 2.59-2.69 vs 2.49-2.62; axial ell3
2.96-3.06 vs 2.75-2.92 (~4-7% higher at at0.1); rho 1.12-1.36 vs 1.04-1.31 (same falls-with-g^2-and-Nf).
at0.1 PNGs added: renorm_dim_fermion_Nf{2,4,6}_at0.1, renorm_rho_at0.1 (shared/renorm_dim now 17).
CAVEAT: L1 at0.1 stats PRELIMINARY (ncfg 456-960) -> errors shrink after 6d topup, re-dump then.
ALIGNMENT: all 324 operator rows config-aligned correlated jk (aligned=Y). rho 29/36 aligned; 7 fall to
independent errors (L3 Nf4 g1.5/3/4.5, L3 Nf6 g1.5/3/4.5, L4 Nf6 g6 -- glue coverage ends earlier than
fermion conn there). central_sys-null glue rows RESOLVED (a7 re-dump).
PNGs (13) in shared/renorm_dim/: renorm_dim_fermion_Nf{2,4,6}_at0.2, renorm_dim_glue_Nf{2,4,6}_at{0.2,0.1},
renorm_dim_perm_Nf2_at0.2, renorm_rho_at0.2.
HEADLINE (Delta = 2 m_O/m_J, at0.2): scalar PS ell0 ~1.75-1.96 (BELOW the current; falls w/ g^2 and L);
scalar ell1 ~2.5-2.7, ell2 ~2.8-3.6, ell3 ~2.9-3.6; axial ell2 ~2.5-2.9, ell3 ~2.75-2.95; vector ell1
~2.1-2.4 (near 2), ell2 ~2.9-3.3; glue F^2 0++ ~3.6-4.4 (near 4), F l=2 ~2.9-3.3 (near 3). rho =
m_axial(ell1)/m_F(l1): 1.31 -> 1.04 falling w/ g^2 and Nf at L1, ~0.9-1.27 L2/L3, L4 ~1.0-1.4 (6-12% err);
crosses 1 at large g^2 for Nf6.
PENDING: at0.1 fermion vector (1d at0.1 dump), joint re-dump after 6d conn topup.
sp EXCLUDED from renorm dump (NM 2026-09-17): the axial sp VSH is a special perambulator-based study, NOT
the standard production channel -> stays recorded in this memo (see "Axial sp VSH" section) but does NOT
enter the phase-2 renorm dump. No axial_sp dump, no jackknifer row.
L4 UPDATE (SCC configs rsynced 2026-09-16, L4 only): a7 re-runs+finalizes L4 F/Fl2, a6 may extend F^2/F^4
to L4 (pending NM); jackknifer refreshes L4 glue renorm rows + rho on their re-dumps. CONSEQUENCE
(jackknifer): the new L4 configs extend the GLUE slab range beyond the FERMION (axial) conn k_max (~739,
9 slabs) => L4 rho rows drop to INDEPENDENT errors UNLESS the axial conn is ALSO measured on the new L4
configs + d4 re-dumps L4 axial_tp. Within-glue L4 ratios (Fl2/F, Fsq/F) stay aligned (same config set).
=> DECISION for NM: commission the fermion (axial) conn L4 measurement on the new configs to keep L4 rho
correlated, or accept independent-error L4 rho.

## Axial (qed3-d4) -- Conn A
STATUS: DELIVERED 2026-09-14 (current working choices; re-run pending conn-complete). Units = LATTICE
dimensionless a_t*m (divide by a_t for physical). ell=0 excluded (current tower); ell=1,2,3 only.
Fit = weighted-const per jk sample (uncorrelated), plateau[2,7] bin10; stat + sys(window +1 shift) + comb.

Authoritative numbers markdowns (in analysis_axial/, source of truth):
- `axial_final_params_claude.md` -- per-L Hankel ladders + at-overrides, fit windows, procedure, regen cmds.
- `hankel_ell_spectrum_masses_at0.200000_claude.md` -- at0.2, ell=1,2,3, ALL Nf x L1-4 (108 rows).
- `hankel_ell_spectrum_masses_at0.100000_claude.md` -- at0.1; **FINAL = L1 rows only** (PRELIMINARY,
  stats topping up). L2 at0.1 rows present but NOT part of the final share.

Final PNGs copied to `final/shared/axial/` (15): hankel_ell_spectrum_L{1-4}_Nf{2,4,6}_at0.200000 (12) +
hankel_ell_spectrum_L1_Nf{2,4,6}_at0.100000 (3, PRELIMINARY).

Headline excerpt (ground ell=1 T1, a_t*m, comb err) -- L1 where both a_t exist:
| L | gsq | Nf | at0.2 a_t*m | at0.1 a_t*m (prelim) |
|---|-----|----|-------------|----------------------|
| 1 | 0.5 | 2 | 0.3557(5)   | 0.1815(5) |
| 1 | 1.0 | 2 | 0.3346(10)  | 0.1728(6) |
| 1 | 1.5 | 2 | 0.3176(15)  | 0.1643(11) |
(full per-ell/per-ensemble numbers in the source md above.) NOTE: interpreting the at0.1-vs-at0.2
comparison (a_t-artifact vs correct scaling) is d4/NM's call -- recorded here as raw numbers only, no
conclusion asserted by the organizer.

## Axial per-m / rotational-symmetry (qed3-60) -- NM-approved KEEPER 2026-09-14
Per-m mass Delta_m (a_t*m, dimensionless) for the axial tp ell=3 tower, m=-1,0,+1 kept SEPARATE, vs 1/L^2.
Nf2, at0.2, smallest gsq per L (L2 g1.0, L3 g1.5, L4 g2.0). Block-Hankel single-state, per-L ladders
(L2[0,2,4] L3[0,3] L4[0,2], reb1@4 T0=1, bin10), SHIFTED const-fit windows L2[2,8] L3[2,7] L4[2,6];
d4's const_fit verbatim, sys=+1 shift, comb. L1 EXCLUDED (m=+-1 vanish: Y_{3,+-1} nodes on the 12 icos
vertices).
| L | 1/L^2 | m=-1 | m=0 | m=+1 |
|---|-------|------|-----|------|
| 2 | 0.2500 | 0.6144(89) | 0.5872(54) | 0.6085(83) |
| 3 | 0.1111 | 0.6814(80) | 0.6844(66) | 0.6888(85) |
| 4 | 0.0625 | 0.7166(123)| 0.7341(71) | 0.7218(104)|
m=+-1 track each other; m=0 below at L2, above at L4.
Files (analysis_axial/): hankel_perm_delta_vs_L_at0.200000_claude.{png,md} (KEEPER); per-m plateau-fit
hankel_perm_platfit_effmass_at0.200000_claude.png/_masses_...md; state doc findings_mvariational_claude.md.
PNGs copied to final/shared/axial/ (delta_vs_L + platfit_effmass). CAVEAT (qed3-60): the vs-1/L^2 axis uses
DIFFERENT gsq per L (smallest per L) -> mixes lattice spacing + coupling, NOT a fixed-physics continuum line.
PENDING: re-run after L4 stats improve (NM working with a remote agent on L4) -> tightens L4.

### qed3-60 per-m extended to SCALAR PS (2026-09-14, NM-approved; coordinated with 5b)
Per-m mass Delta_m (a_t*m) for the scalar PS ell=3 tower (icos T2+G), m=-1,0,+1 SEPARATE, vs 1/L^2.
Nf2, at0.2, smallest gsq per L (L2 g1.0, L3 g1.5, L4 g2.0). Block-Hankel single-state; scalar PS conn
C=-(2 Re Vpp) (LOOP_SIGN=-1, 5b-verified). L1 excluded (m=+-1 icos-vertex nodes). LOCKED SCALAR config:
off L2[0,3] L3[0,3] L4[0,2]; reb1@4 T0=1; bin10; windows L2[2,6] L3[2,5] L4[2,5] (scalar-tuned, narrower
than axial; T0-stable).
| L | 1/L^2 | m=-1 | m=0 | m=+1 |
|---|-------|------|-----|------|
| 2 | 0.2500 | 0.745(22) | 0.720(11) | 0.736(14) |
| 3 | 0.1111 | 0.806(18) | 0.796(25) | 0.806(21) |
| 4 | 0.0625 | 0.861(27) | 0.843(14) | 0.845(23) |
Per-m tightly clustered per L; mass rises toward smaller 1/L^2. Files: source in analysis_axial/
(qed3-60's per-m tooling) hankel_perm_delta_vs_L_scalar_ps_at0.200000_claude.{png,md} + platfit; the 2
PNGs copied to final/shared/scalar/. Same gsq-per-L axis caveat (mixes spacing+coupling, not fixed-physics).

## Axial sp VSH (qed3-4e) -- 2026-09-17
Tangential sp = VSH decomposition (electric Phi_l, magnetic Psi_l) vs tp radial, Nf2 gsq1.0 at0.2,
dimensionless a_t*m, GPOF reb1@4 T0=1, from COMPLETE(L1)/truncated-smearing-benign(L2) distillation perams.
Truth = final/analysis_axial/findings_sp_vsh_claude.md sec 7. tp reproduces d4's tp tower (L2 ell1,2,3 =
0.3664/0.5075/0.6001) -> pipeline VALIDATED; V=A everywhere.
- L2 (627 cfg): ell1 Phi_1 0.3679(5) Psi_1 0.4978(48) tp 0.3664(7); ell2 Phi_2 0.4995(86) Psi_2 0.3695(16)
  tp 0.5072(20); ell3 Phi_3 0.3706(99) Psi_3 0.4186(155) tp 0.5967(122).
- L1 (398 cfg): ell1 Phi_1 0.3321(18) Psi_1 0.3604(41) tp 0.3332(14); ell2 Phi_2 0.3471(47) Psi_2 0.3327(28)
  tp 0.4175(17).
PHYSICS (qed3-4e): (1) the low higher-ell tangential values are NOT aliasing (VSH Gram overlaps orthogonal
to 1e-16, free continuum gives clean Phi=Delta+l-1/Psi=Delta+l) -- they are a LARGE a^2 artifact SPECIFIC
to the tangential (curl/grad on the simplicial sphere): dim<=3 all channels agree, but at dim4 radial tp3
~0.60 (near-continuum) while tangential Phi_3,Psi_2 ~0.37 collapse; heals L1->L2 => fix is finer L / a->0,
NOT a fit change. (2) degeneracy Psi_l ~ Phi_{l+1} (same continuum dim Delta+l AND same lattice value:
dim3 0.498 vs 0.500, dim4 0.370 vs 0.371). tau_gw adjoint-leg bug (flagged by Fin:Two-meson): sp axial legs
UNAFFECTED (conj-transpose-of-forward AblkS = delta-tau, correct).

## AXIAL at0.1 Hankel param reference (for relaying to scalar/vector when asked)
at0.1 DIFFERS from at0.2 via d4's ATOFF/ATWIN overrides (axial_final_params_claude.md):
- L1: at0.2 offsets [0,3,6] win [8,15]  ->  at0.1 offsets [0,4,8] win [12,24] (flatter ladder + long
  clean plateau at at0.1).
- L2: [0,2,4] win [4,10] -- SAME at both a_t.
- L3/L4: no at0.1 data.
Common both a_t: reb1@4, T0=1, kmin=20, bin10, dimensionless a_t*m. (Per-channel FINAL at0.1 params are
still NM's to set; this is the axial reference to mimic.)

## Vector (qed3-1d)
STATUS: IN PROGRESS 2026-09-14 -- NM released vector; 1d pinged, organizer directed her. Mimic the
HANKEL ell-spectrum path (copy d4's latest Hankel core+driver, minimal edits; 5b's scalar adaptation is
the worked template). Loader = vector physical correlator with disc DC+plateau subtraction folded in
(full = conn + (sub - plat); conn = -gl_conn m-avg no factor 2). ell=0 DROPPED (vector current tower) ->
ell=1,2,3. Same fit (weighted-const per-jk plateau, sys=window+1, comb), dimensionless a_t*m, m-avg,
schema L|gsq|Nf|ell|irrep|off|win|a_t*m|stat|sys|comb. GATE: NM supplies vector Hankel params
(OFFMAP/REBT/NKEEP/T0) + fit windows (WINMAP); 1d writes impl_plan, builds pipeline, then asks NM.
On final delivery: record md here + copy PNGs to final/shared/vector/.
- 2026-09-14 UPDATE: pipeline BUILT + VALIDATED (analysis_vector/): effmass_vector_hankel_core_claude.py
  (shared core verbatim + load_vector folding in disc DC[Eq.12]+plateau[16,Nt/2] -> physical -C_conn+disc_sub)
  + hankel_ell_spectrum_vector_claude.py (5b's driver, ELLS=[1,2,3] ell=0 dropped) +
  vector_hankel_spectrum_impl_plan_claude.md. Validated L1 Nf2 g1.0 (metric PD; ell=1 small-t ground
  a_t*m ~0.32-0.35 matches old mislabeled physical 1.776x0.2=0.355). KEY finding for NM's window choice:
  vector signal clean ONLY at SMALL t (t<=4), disc noise kills t>~5 => vector fit window is small-t
  (~[2,4]/[2,5]), NOT axial's [8,15]. At NM gate (asking OFFMAP/ATOFF/REBT/NKEEP/T0MAP + WINMAP/ATWIN).
- 2026-09-14 DELIVERED (PROVISIONAL -- NM: "same fit params for now, retune later"; used axial params):
  Numbers md (analysis_vector/): hankel_ell_spectrum_vector_masses_at0.200000_claude.md (108 rows L1-4)
  + _at0.100000_claude.md (36 rows: L1 all-g + L2 g2.0). cols L|gsq|Nf|ell|irrep|off|win|a_t*m|stat|sys|comb.
  PNGs copied to final/shared/vector/ (18): hankel_ell_spectrum_vector_L{1-4}_Nf{2,4,6}_at0.200000 (12)
  + _L{1,2}_Nf{2,4,6}_at0.100000 (6). READ (ell=1 T1): L2 (win [4,10]) CLEAN, a_t*m ~0.28-0.46 (e.g.
  L2 Nf2 g1.0 = 0.348(21) -> phys 1.74 ~ old 1.78); L1 (win [8,15]) is NOISE -- axial window too late,
  vector L1 needs small-t ~[2,4] (RETUNE PENDING NM). L3/L4 disc-limited (<100 cfg) provisional; reliable
  full conn+disc = L1+L2. STATUS: PROVISIONAL pending NM window retune (esp L1).
- 2026-09-14 STATE-SAVE (1d, diagnostics; production driver NOT yet repointed -> no final md/PNGs in
  shared/vector/ yet). LOCKED by NM ("use for now"): Dt=[0,2], T0=2, reb1@3, NKEEP=1, bin10, kmin20.
  ELLS=[1,2] (ell=0 conserved-current dropped; ell=3 DROPPED by NM). Windows L1 [5,10], L2 [4,8].
  Reliable = L1+L2 (L3/L4 disc-limited). Method = block-Hankel/GPOF on combined -C_conn + (disc DC[Eq.12]
  + plateau[16,64] sub), added per-config, m-avg ground, dimensionless a_t*m.
  RESULTS (a_t*m at0.2, comb): ell=1 T1 L1 ~0.395-0.403, L2 ~0.369-0.403 (comb ~0.006-0.025); ell=2 H
  L1 ~0.49-0.53, L2 ~0.48-0.54 (comb ~0.014-0.042). ell=1 ~Nf-independent but NM SET ASIDE Nf-dep here.
  Cross-checks: matches old (mislabeled-physical) vec mass + plain cosh effmass; mild L1(~0.40)>L2(~0.38)
  = finite-vol/spacing. KEY: the double-subtracted disc piece is CLEAN (not the noise source); the earlier
  "terrible" masses were a fit-WINDOW artifact (axial [8,15] sits in vector noise). Deliverable so far:
  ell1_vs_invNf_vector_at0.200000_claude.{png,md} (in analysis_vector/, NOT yet copied to shared). 1d will
  repoint the driver to Dt=[0,2]/these windows + re-emit final md/PNGs when NM greenlights -> then copy.
- 2026-09-15 VECTOR at0.2 FINAL (1d; driver repointed to locked params). Numbers md refreshed:
  analysis_vector/hankel_ell_spectrum_vector_masses_at0.200000_claude.md. 6 refreshed L1,L2 PNGs re-copied
  to final/shared/vector/ (L{1,2}_Nf{2,4,6}_at0.200000). Locked: Dt=[0,2] T0=2 reb1@3, ELLS=[1,2], windows
  L1[5,10]/L2[4,8], LEXICOGRAPHIC sorted(glob) binning (matches d4/5b fermion). FINAL numbers a_t*m at0.2:
  ell1(T1) L1~0.40 L2~0.38; ell2(H) L1~0.51-0.53 L2~0.48-0.53. Reliable = L1+L2 (physical conn+disc).
  shared/vector/ cleanup (1d decision 2026-09-15): KEEP the 6 at0.1 PNGs (pre-lock; 1d re-emits at0.1 with
  locked params when that arm runs -> overwrites in place; marked PENDING). PRUNE the 6 L3,L4 at0.2 PNGs
  (pre-lock, disc-limited, out of reliable scope, NOT regenerated -> permanent stale) -- flagged to NM to
  run the rm (organizer does not delete). FINAL vector full conn+disc = L1+L2 only. Phase-2 jk dump: shared/jk_dumps/vector/vector_at0.2_claude.json (validated).
  RENORM ALIGNMENT CAVEAT (1d->jackknifer+NM): vector config set = conn INTERSECT disc (fewer configs than
  axial_tp conn-only reference) -> bin_k won't exactly match axial_tp; needs a disc-matched reference subset
  or k-interval re-bin (jackknifer's call; independent-error fallback otherwise).

## Scalar (qed3-5b)
STATUS: IN PROGRESS 2026-09-14 -- directed by NM to MIMIC d4's axial procedure (per-L Hankel ell-spectrum,
dimensionless a_t*m, kmin=20, weighted-const per-jk fit, stat+sys(window+1)+comb, bin10) and generate the
equivalent scalar figures. SCALAR KEEPS ell=0. NM will supply the DETAILS (scalar fit ranges + Hankel
parameters); 5b sets up the pipeline then asks NM directly for those. On final delivery: record numbers md
here + copy PNGs to final/shared/scalar/.
- 2026-09-14 UPDATE: pipeline BUILT (copied d4's latest scripts, minimal edits) -- core
  effmass_scalar_hankel_core_claude.py (verbatim Hankel; loader -> PS=2ReVpp, FS=Re(Vpp+Vmm^FS)) +
  driver hankel_ell_spectrum_scalar_claude.py (ELLS=[0,1,2,3] keeps ell=0; --chan PS|FS|both), in
  analysis_scalar/. Smoke-tested OK (ell=0 PS==FS ~0.317; ell=1 ~0.42 matches old physical PS). NOT
  deliverable yet (placeholder params). AT ask-NM gate: requested Hankel params (OFFMAP/REBT/NKEEP/T0),
  fit windows (WINMAP), + loop-sign confirm. Plan: at0.2 (L1-4) + at0.1 (L1,L2) PS+FS.
- PHYSICS FLAG (5b -> NM): stored 2 Re Vpp is NEGATIVE for t>=1 (missing closed-fermion-loop sign) ->
  Hankel metric non-positive; physical C = -(2 Re Vpp) is reflection-positive + reproduces old PS mass,
  so 5b applies LOOP_SIGN=-1 (verified, flagged to NM for confirm).
- 2026-09-14 SCALAR PS at0.2 FINAL (NM fixed the params). PS ONLY -- confirmed PS==FS (GW no-op; direct
  per-config check median ~1e-5), FS dropped. ell KEPT 0,1,2,3 (0=A 0++ ground, 1=T1, 2=H, 3=T2+G).
  Method = d4 block-Hankel core VERBATIM (loader swapped); reb1@4 T0=1; offsets L1 0-3-6 / L2 0-2-4 /
  L3 0-3 / L4 0-2; dimensionless a_t*m; kmin=20 bin10; weighted-const per-jk + sys(+1) + comb.
  FINAL windows (scalar plateaus onset ~t3): L1 [3,10], L2 [3,9], L3 [5,10], L4 [3,8]. LOOP_SIGN=-1.
  Numbers md: analysis_scalar/hankel_ell_spectrum_scalar_PS_masses_at0.200000_claude.md (144 rows,
  L1-4 x Nf{2,4,6} x ell{0,1,2,3}). Params: analysis_scalar/scalar_final_params_claude.md.
  PNGs copied to final/shared/scalar/ (12): hankel_ell_spectrum_scalar_PS_L{1-4}_Nf{2,4,6}_at0.200000.
  QUALITY: ell=0 (0++ ground) + ell=1 EXCELLENT all L; ell=2 good; ell=3 clean L1, LOW-CONF L2/L3
  (excited-state bump, large sys). ell=0 a_t*m ~0.29-0.37 (falls with g^2, rises with L).
  2026-09-15 UPDATES (5b, from Phase-2 dump prep): (a) baked the FINAL windows [3,10]/[3,9]/[5,10]/[3,8]
  into the driver WINMAP (had been CLI overrides; a stale module default [8,15] was briefly read then
  fixed -- published table was always the correct [3,10], re-verified dump==table 0 mismatch). (b) L3 Nf2
  g4.5 GREW 373->430 cfg (6d topup) -> recomputed that ensemble fresh at 430; refreshed the at0.2 table
  (that ensemble's 4 rows) + the L3_Nf2 at0.2 PNG (re-copied to shared/scalar/); other 11 PNGs + rows
  unchanged. NOTE: scalar+axial should RE-DUMP together after the 6d topup settles so bin_k align on
  grown ensembles (transient mismatch -> jackknifer independent-error fallback).
  PENDING (not blocking): possible dedicated ell=3 window at L3.
- 2026-09-14 SCALAR PS at0.1 L1 arm DONE (PRELIMINARY -- L1 stats topping up; L2 DEFERRED by NM).
  Params: canonical at0.1 convention -- L1 offsets [0,4,8], reb1@4, T0=1 (= d4 axial); const-fit window
  [12,24] KEPT AS-IS (NM: tough call, keep for now; NOT scalar-retuned -- scalar at0.1 plateaus long+flat
  from t~4 so an earlier onset e.g. [6,24] is available -> possible later revisit). Numbers md:
  analysis_scalar/hankel_ell_spectrum_scalar_PS_masses_at0.100000_claude.md (36 rows, L1 x Nf{2,4,6} x
  ell0-3). PNGs copied to final/shared/scalar/ (3): hankel_ell_spectrum_scalar_PS_L1_Nf{2,4,6}_at0.100000.
  a_t-ARTIFACT CHECK (5b's read, recorded as raw obs): ell=0 a_t*m ~0.16 at at0.1 vs ~0.32 at at0.2 =>
  same physical m ~1.6 (dimensionless halves, physical ~constant). Masses-md headers now at-AWARE (print
  per-L off/window per at); at0.2 rows byte-identical (no PNG re-copy needed).

## Glue -- SPLIT: qed3-a7 (linear F l=1 + Fl2 l=2) / qed3-a6 (F^2/F^4 0++)
STATUS: at0.2 baseline complete (L1-L4 all OK). a7 F l=1 at0.2 FINAL + a6 F^2/F^4 0++ at0.2 FINAL (below).
a7 next when NM directs: at0.1 F, then Fl2 (l=2).

- 2026-09-14 GLUE F^2/F^4 0++ (qed3-a6 "Fin: Fsq") at0.2 FINAL (NM signed off). Numbers md:
  analysis_gluonic/Fsq_final_masses_at02_claude.md. PNGs copied to final/shared/glue/ (3):
  Fsq_spectrum_at02, Fsq_grid_L1_at02, Fsq_grid_L2_at02. METHOD (locked): shape-basis GEVP, 14 ops
  (7 icos shapes x p2[F^2]+p4[F^4], face_sign ON, l=0) + EXPLICIT identity op (add_const=1) + vacsub=1
  => state0 = exact vacuum, 0++ = STATE 1. Fixed-t0 t0=0 tr=1, NO Hankel (rejected: identity row
  duplicated across offsets -> singular). Dimensionless a_t*m (AT=1.0), kmin=20, binsize=80, window [1,6]
  both L, diagonal-const per-jk + sys(+1) + comb. Own binary glue_gevp_analysis_fsq_claude.o (add_const
  arg27, default-off = bit-identical to a7's glue_gevp_analysis_t0_claude.o; a7's binary UNTOUCHED).
  RESULTS (a_t*m 0++): L1 ~0.51-0.56 all Nf/gsq (near free anchor 2sqrt2*0.2=0.566, ~Nf-flat; L1
  chi2/dof 2-7 as effmass droops past t~4). L2 RISES with gsq AND Nf: Nf2 0.58->0.66, Nf6 0.62->0.75
  (gsq 1->3), chi2/dof 0.6-2.0. Scope = L1+L2 x Nf{2,4,6} x per-L gsq (18 ens), at0.2 massless.
  Boundary held (F^2/F^4 0++ only; no Fl2, no sigma^2 mixing).
  L4 RESOLUTION 2026-09-17 (ToFin): 799 IS the FULL L4 set (target bumped 600->800, 799=target-1=DONE;
  SCC stats already landed weeks ago, Nf4/6 640-756 -> 799). So a7's 799-config L4 F/Fl2 ARE full-stats
  final (a7 to CONFIRM analysis read current 799). Gauge configs = SIBLING dir Nf<n>_gsq<g>at0.2...L4_hb.../
  ckpoint_lat.<k> (not data_). ONE gap: Nf6 g6.0 = 778/799 glue h5 (21 unmeasured) -> ToFin preps CPU
  topup glue_msm_shapes.
  CLOSED 2026-09-17 (a7 confirmed): L4 F l=1/Fl2 is FULLY FINAL for all 9 ensembles as-is -- the 8
  non-Nf6g6 read the full 799 (Nc=780, 9 bins). Nf6 g6 topup (778->799) is a NUMERICAL NO-OP: the 21
  missing configs are the TAIL k=779-799, which fall in the DROPPED binsize-80 jackknife remainder (fit uses
  first 720 = k=20..~739), so nbins stays 9 and the number is identical before/after. => NO re-run, NO
  re-dump, NO measurement needed; all L4 F/Fl2 numbers + glue_F/glue_Fl2 jk-dumps STAND as final. Topup is
  OPTIONAL h5-completeness hygiene only (skipped by default). (To actually add a bin at Nf6 g6 would need
  >=800 post-kmin = >819 ckpts, above the 800 ceiling.) OPEN:
  if SCC generated beyond the 800 target, a fresh rsync_scc / SCC-agent maxk would confirm (flagged to NM);
  else 799 = ceiling.
  2026-09-16: a6 ran L3 (999 cfg) + L4 (799 cfg) 0++ through the locked pipeline as a diagnostic -> NO
  clean plateau (0++ effmass droops from t=1 into noise by t~2-3; strongest coupling pure noise). NOT
  nops-limited (14-op l=0 basis at every L) -- intrinsic 0++ S/N at coarse temporal resolution. NM
  DECISION: F^2/F^4 0++ STAYS L1+L2 ONLY. No L3/L4 md/PNG refresh, no L4 jk-dump; glue_Fsq at0.2/at0.1
  dumps stand unchanged under v1.1. (=> the SCC L4 update affects ONLY a7's F l=1 + Fl2, not a6's 0++.)
- 2026-09-14 GLUE F^2/F^4 0++ at0.1 FINAL (qed3-a6; NM signed off). Numbers md:
  analysis_gluonic/Fsq_final_masses_at01_claude.md. PNGs copied to final/shared/glue/ (3):
  Fsq_spectrum_at01, Fsq_grid_L1_at01, Fsq_grid_L2_at01. Same locked method as at0.2 (14-op shape GEVP +
  identity + vacsub=1 -> 0++=state1, fixed-t0 t0=0 tr=1 no-Hankel, AT=1.0, binsize=80, kmin=20); flow-at
  does NOT gate masses so at0.1 h5 usable as-is (no re-dump). nbins=24. Windows L1[1,7] L2[1,5] (at0.1
  plateau longer than at0.2), chi2/dof 0.08-2.0. Coverage = 12 ens (L1 Nf{2,4,6} g{0.5,1,1.5} + L2
  Nf{2,4,6} g2.0 = all that exists at0.1). RESULTS (a_t*m 0++): L1 ~0.255-0.30 (on at0.1 anchor 0.283,
  ~Nf-flat); L2 g2.0 rises w/ Nf (Nf2 0.318 / Nf4 0.335 / Nf6 0.351). a_t*m(at0.1) ~ HALF a_t*m(at0.2)
  => consistent physical m ~2.7-2.8 at both a_t (mild a_t artifact).

- 2026-09-14 GLUE F l=1 (linear F_12, T1 glueball) at0.2 FINAL. Numbers md:
  analysis_gluonic/F_final_masses_at02_claude.md (all 36: L1-4 x Nf{2,4,6} x gsq). PNGs copied to
  final/shared/glue/ (7): F_fixt0_spectrum_at02 (KEY a_t*m vs g^2 per L, Nf2 red o/Nf4 blue s/Nf6 green ^),
  F_fixt0_grid_L{1,2,3,4}_at02 (per-L effmass montages, m-avg + fitted plateau), F_binsize_scan_Nf2g0.5L1
  + F_binsize_scan_multi. Method (consistent w/ other channels): FIXED-t0 GEVP t0=0 tr=2, NO Hankel (Hankel-on-
  shape-matrix shelved -- inflates/singular for glue), m-avg over l=1 m-triplet, kmin=20, dimensionless
  a_t*m (phys = /0.2), BINSIZE=80 (autocorr-saturated). Windows L1[1,12] L2[1,8] L3[1,5] L4[1,5].
  Impl: gevp_t0_variants_impl_plan_claude.md; binary glue_gevp_analysis_t0_claude.o (t0mode/t0fix/offs;
  sliding default bit-identical).
  Representative a_t*m(err): L1 Nf2 g0.5 0.2717(31); L2 Nf6 g3 0.3709(93); L3 Nf6 g4.5 0.4315(81);
  L4 Nf2 g2 0.3061(134). TREND (NM key): within each L, F RISES with Nf (Nf6>Nf4>Nf2) AND with gsq;
  chi2/dof mostly <2. CAVEAT: L4 (~9-10 jk bins, Nc~780) binsize does NOT saturate -> L4 errors are
  ~LOWER BOUNDS (stats-limited); L1/L2/L3 clean.
- 2026-09-14 GLUE F l=1 at0.1 FINAL (a7; NM LOCKED windows L1[1,12] L2[1,8]). All 12 at0.1 ens (L1
  g{0.5,1,1.5} + L2 g2.0, all Nf) -- measurement GAP FILLED. Same recipe as at0.2 (fixed-t0 t0=0 tr=2,
  no Hankel, m-avg, kmin=20, binsize=80, dimensionless a_t*m). Numbers md:
  analysis_gluonic/F_masses_at01_claude.md (header FINAL). PNGs copied to final/shared/glue/ (4):
  F_at_check_phys (KEY physical-m a_t comparison), F_fixt0_spectrum_at01, F_fixt0_grid_L1_at01,
  F_fixt0_grid_L2_at01. RESULT (a7's read): physical F l=1 ~a_t-INDEPENDENT (at0.1 ~ at0.2, mild few-%
  down at at0.1) -> glue F scales properly, UNLIKE the fermion a_t-artifact; a_t*m at0.1 ~0.13-0.16
  (~half at0.2). => BOTH glue F l=1 at0.2 AND at0.1 now FINAL.
- 2026-09-15 GLUE Fl2 (l=2, H irrep) at0.2 FINAL (a7; NM: ell=2, no at0.1 for l=2). Numbers md:
  analysis_gluonic/Fl2_masses_at02_claude.md (36 ens). Method = same as F l=1 (fixed-t0 t0=0 tr=2, no
  Hankel, m-avg over 5-fold l=2 H multiplet, kmin=20, binsize=80, dimensionless a_t*m); LOCKED windows
  L1[1,5] L2[1,4] L3[1,4] L4[1,3]. PNGs copied to final/shared/glue/ (6): Fl2_fixt0_spectrum_at02,
  Fl2_over_Fl1_ratio_at02, Fl2_fixt0_grid_L{1,2,3,4}_at02. RESULT: a_t*m ~0.43 (L1) -> ~0.61 (L3 Nf6
  g4.5); same Nf+gsq trend as Fl1; Fl2/Fl1 ~1.5-1.6. Caveat: L1 gsq0.5 chi2/dof=4.99 (only poor fit;
  [1,4] would clean it). => GLUE F sector FINAL: F l=1 (at0.2+at0.1) + Fl2 l=2 (at0.2). ell=3 DEFERRED
  by NM (needs re-dump + T2/G irrep split).
