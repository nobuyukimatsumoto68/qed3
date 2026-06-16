# Diluted exact-K + local current-current correlator -- impl plan

File: **`jj_corr_dilute_claude.cu`** (verbatim copy of `jj_corr_block_t_claude.cu`; dilution + local
currents + origin superposition added per the chunks below). Companion: `jj_local_stoch_estimator_design_claude.md`
(estimator algebra), `jj_stoch_local_impl_plan_claude.md` (local conventions), `project_state.md` (2026-06-12).

**Algorithm sources (cite in code header):**
- dilution / all-to-all propagators -- Foley, Juge, O'Cais, Peardon, Ryan, Skullerud, "Practical all-to-all
  propagators for lattice QCD," hep-lat/0505023.
- source-origin superposition / translation (volume) averaging = the **master-field** approach -- orig.
  M. Luscher, "Stochastic locality and master-field simulations of very large lattices," arXiv:1707.09758
  (EPJ Web Conf. 175 (2018) 01002); algo follow-up Francis, Fritzsch, Luscher, Rago, arXiv:1911.04533;
  framing used by us in arXiv:2301.08696.

---

## 1. Goal & scope
ONE program computing, in a single pass that **shares the diluted `phi'=D_m^{-1}eta`**:
- **exact conserved current K** -- connected **tp + sp** + **disc** + **axial** (jj_corr_block_t content, minus ylm);
- **local (ultralocal) current** -- connected **sp/tp** = diagonal `s1,s2,s3` + **local disc** `J_a(t)`.
With **spin x even/odd-time dilution = 4 patterns** AND **source-origin superposition** (Sec. 3.5).

**Why merged (supersedes "separate local from K"):** `phi'=D_m^{-1}eta` is the COMMON sink leg; exact-K and
local differ only in the *insertion* (`K` vs point-`sigma`). The local sink is a cheap local `sigma`-multiply
on the shared `phi'` (no extra solve); local adds only its own source legs. One pass emits both, sharing the
diluted `phi'` + eta + 4-pattern machinery + output.

**OUT of this file:**
- **ylm** -- LOCAL-current observable in a SEPARATE one-end-trick file (`jj_local_ylm_stoch`, ~copy of
  `meson_pq_wall_v2`; the $Y_{lm}$ double-sum is a different method).
- **local axial** -- future TODO (GW $(1-D_{ov})$ care; vector first). Exact-K axial IS included.

**Validation (free, L=1):** exact tp/sp vs `corr_deter_exactsum_L1`; local s1,s2,s3 vs `corr_deter_local_L1`
(both deterministic, zero-noise), reading off the dilution variance reduction. Production
`jj_corr_block_t_claude.cu` untouched (A/B).

---

## 2. The dilution scheme (4 patterns)
Patterns `p=(t_s,spin)`, `t_s in {0,1}`, `spin in {0,1}`:
```
eta.time_spin_dilution(rng, t_s, /*t_block=*/Comp::Nt/2, spin);   // interval = Nt/t_block = 2
```
- `t_block=Nt/2=64` => `interval=2`: `t_s=0`->even, `t_s=1`->odd; `spin` fixes the component.
- 4 DISJOINT supports tile the volume => diluted estimator = plain sum over the 4 patterns (unbiased).
- Removes cross-class noise: dominant `t0<->t` cross-term for OPPOSITE-parity separations (~half) + spin
  off-diagonal cross-terms (spin dilution also supplies the $\sigma_1,\sigma_2$ cross-spin terms, both exact
  and local, as in `disc_claude.cu`'s `accumulate_loop_gamma`).
- Caveat: `K(a,t0)` / local temporal hop span `{t0,t0+1}` (opposite parities) => STANDARD 4-class sum,
  partial reduction; NOT the "only-t0's-class" free reduction. Finer time slicing = next rung.

---

## 3. Estimator (what shares what)
Per hit, per pattern `p`, the SHARED solve `phi'^(p)=D_m^{-1}eta`, then:
```
# exact K (summed-origin source, ONE solve per insertion -- Sec. 3.5):
psi_tp[n]^(p) = D_m^{-dag} ( sum_b K^dag(n ,t0_b) eta )        # sites
psi_sp[a]^(p) = D_m^{-dag} ( sum_b K^dag(lk,t0_b) eta )        # links
exact sink: K(ins,t) phi'^(p)  (block_t)  -> exact conn tp/sp ; exact disc rides (eta^dag K(ins,t)phi')
# local sp/tp (summed-origin source, ONE solve per (channel,site)):
chi[c,n]^(p) = D_m^{-dag} ( sum_b sigma_c(n,t0_b) eta )        # c in {1,2,3}, sites n
local sink: sigphi_c = sigma_c phi'^(p)  (LOCAL multiply, no solve)
  conn_c(t) += w_site[n] * chi[c,n]^(p)(t,n)^dag . sigphi_c(t,n)     # SUPERPOSED over origins
  local disc J_c(t) += w_site[n] * eta(t,n)^dag sigphi_c(t,n)        # rides phi'
```
- E[sum_p ...] reproduces the (superposed) diagonal traces with lower variance; signs/normalization match
  the deterministic refs (exact->`corr_deter_exactsum`, local->`corr_deter_local`).
- **Weights:** exact `w_tp=dual_areas/kappa_t^2`, `w_sp=link_volume/kappa^2`; local `w_site=dual_areas`
  (no 1/kappa^2). `1/(4pi)` folded by write_corr; discs RAW.

### 3.5 Source-origin SUPERPOSITION = master-field (replaces the t0 loop)
This is the **master-field** approach (translation/volume averaging from few sources; orig. Luscher
arXiv:1707.09758; our framing in arXiv:2301.08696). Instead of looping over `n_t0` origins (separate source
solves + per-t0 output), insert at BOTH origins `t0 in {0, Nt/2}` SIMULTANEOUSLY: the source is their SUM
(above), one solve per insertion. The measured correlator is the superposition
```
C(t) = G(t) + G(t - Nt/2)   (no cross-term in the mean; periodic, period Nt/2; symmetric about Nt/4)
```
**Analysis fold:** use `C(t)` on `[0, Nt/2)`, fold about `Nt/4`; the far tail `G(t+Nt/2)` is the (small,
fast-decaying) WRAP-AROUND contamination -- "a lot of wrap-around, but fine for our purpose" (signal dies
well before Nt/4). Store the origin set `t0s={0,Nt/2}` in the .h5 for the fold. (`n_t0` becomes "number of
superposed origins"; default 2. Source legs NO LONGER scale with `n_t0`.)

---

## 4. Cost & validation
- **Per pattern:** `phi'` (1) + exact source `n_sites+n_links` (ONE summed-origin solve each) + local source
  `3*n_sites` (ONE each) + sinks (exact block_t; local = free multiply). **x4 patterns.** No `n_t0` factor
  on the source legs now (superposed). At L=1: `4*(1 + 42 + 36) = 316` solves/hit.
- **Fair test:** 1 diluted hit vs 4 plain hits (equal cost) -- expect lower variance.
- **Validation notebook:** copy/extend `jj_free_stoch_validate_claude.ipynb` -- exact tp/sp on
  `corr_deter_exactsum_L1`, local s1,s2,s3 on `corr_deter_local_L1`, with diluted-vs-volume error bars at
  matched cost, AND the fold of the superposed `C(t)`.

---

## 5. Output naming (per HIT; superposition removes the t0 split)
`data_<ESNID>/corr_dil_nt0<N>_nhits<H>/corr.<k>.h<h>.h5` (one file per config k, hit h). The earlier
per-`(h,t0)` split is SUPERSEDED -- origins are superposed in the source, so there is no separate `t0`.
- Keys: exact `h0/{tp,sp}/{Vpp,Vmm}`, axial `Apm`, `h0/disc/{tp,sp}/J`; local `h0/{s1,s2,s3}/{Vpp,Vmm}`,
  `h0/disc/{s1,s2,s3}/J`; plus `t0s={0,Nt/2}` (for the fold) and `rng_seed`. (Drop the `t0_<b>` nesting.)
- Resume = skip a completed `(h)` file (own `complete`; `.tmp`+rename per h).
- **Post-dilution default `--nhits 1`** => a single `corr.0.h0.h5` (2 origins folded in).

---

## 6. RNG seeding from string (stable, stored in the .h5)
Canonical = `std::seed_seq` (deterministic & portable; NOT `std::hash`). No `rng.h` change.
```cpp
#include <random>
#include <cstdint>
static int seed_from_string(const std::string& s){       // file-scope helper
  std::seed_seq seq(s.begin(), s.end());
  std::uint32_t w; seq.generate(&w, &w + 1);
  return static_cast<int>(w);
}
...
const std::string seed_str = esnid + "_k"+std::to_string(k) + "_h"+std::to_string(h);   // per (config,hit)
Rng rng(base, seed_from_string(seed_str));
h5.createDataSet("rng_seed", seed_str);
```
Seed per `(k,h)` (the t0 subtlety is moot now -- origins are superposed in one `eta`). Reproduce via
`std::string s; f.getDataSet("rng_seed").read(s); Rng(base, seed_from_string(s));`.

---

## 7. Implementation chunks (line refs into the copy)
- **Chunk 0 (DONE):** verbatim copy `jj_corr_dilute_claude.cu`.
- **Chunk 1 -- header + output tag + seed helper.** Header comment: exact+local diluted+superposed variant
  + cite hep-lat/0505023; `dir_out` "corr_nt0" -> "corr_dil_nt0" (~360); add `seed_from_string` +
  `#include <random>,<cstdint>`. Files: jj_corr_dilute_claude.cu.
- **Chunk 1b -- strip ylm.** Comment out (don't delete): `N_ELL`/`W_ell` (~351), `psi_yl`(375)+`IYL`(383),
  ylm source solve (495), ylm sink + output. Files: jj_corr_dilute_claude.cu.
- **Chunk 2 -- accumulators (single, superposed).** Collapse the per-`b` correlator arrays to ONE each
  (superposed): exact Ctp/Csp, axial Apm, disc Jtp/Jsp; local (NEW) Cs1/Cs2/Cs3 + disc Js1/Js2/Js3. Zero
  per hit; sum over the 4 patterns. Files: jj_corr_dilute_claude.cu (~437-).
- **Chunk 3 -- superposed source + pattern loop + seed.** Reseed per hit
  (`Rng(base, seed_from_string(esnid+"_k"+k+"_h"+h))`, replacing fixed `Rng(base,1234)` line 260). Replace
  `eta.fill_z2_source` (438) with `for(t_s in{0,1}) for(spin in{0,1})`: `time_spin_dilution`; `phi'`
  (462-463). Build each exact source as the SUM over origins `sum_b K^dag(ins,t0_b) eta` (one solve;
  collapse the `for b` loop 477- into a summed RHS), sink (500-), accumulate the single superposed buffers.
  Files: jj_corr_dilute_claude.cu.
- **Chunk 3L -- local current (same pattern loop).** Add `w_site[n]=dual_areas` (near 332-338). Per pattern,
  per (channel c, site n): summed-origin localized source `sum_b sigma_c(n,t0_b) eta` (restrict eta to site
  n at each t0_b, `mult_sigma(c)`, sum), solve `chi[c,n]=D_m^{-dag}(...)`; local sink `sigphi_c=sigma_c phi'`;
  accumulate `Cs{c} += w_site[n]*chi(t,n)^dag sigphi_c(t,n)` + disc `Js{c}[t] += w_site[n]*eta(t,n)^dag sigphi_c(t,n)`.
  (Conventions: `jj_local_deter_claude.cu` + `jj_local_stoch_estimator_design_claude.md`.) Files: jj_corr_dilute_claude.cu.
- **Chunk 4 -- write per hit + t0s + rng_seed.** Fold 1/4pi + write `corr.<k>.h<h>.h5`: exact (tp/sp/disc/axial)
  + local (s1,s2,s3/disc) keys (no t0 nesting) + `t0s={0,Nt/2}` + `rng_seed` + `complete`; resume + `.tmp`-rename
  per h. Files: jj_corr_dilute_claude.cu.
- **Chunk 5 -- run + validate.** `tmp_*_claude.sh` (free L=1, a few diluted hits) + validation notebook
  (exact vs exactsum, local vs corr_deter_local, the superposed-`C(t)` fold, diluted-vs-volume). Files: tmp, notebook.

---

## 8. Resolved decisions (2026-06-12)
1. **Merge local into this file** -- YES (shared `phi'`); supersedes "separate local from K".
2. **Disc** -- exact + local disc, SAME 4 dilution patterns (ride `phi'`).
3. **Axial / mass / parity** -- exact-K axial kept + diluted; mass/parity kept. Validation focus = vector
   (exact tp/sp + local s1,s2,s3). Local axial = future TODO.
4. **ylm** -- OUT (separate local one-end-trick file).
5. **Seeding** -- string -> seed_seq -> int (per (k,h)); `rng_seed` stored per file.
6. **Source-origin SUPERPOSITION** -- insert at {0, Nt/2} simultaneously (summed source, one solve per
   insertion); superposed `C(t)=G(t)+G(t-Nt/2)`; analysis folds; wrap-around accepted. REMOVES the t0 loop
   and the per-`(h,t0)` output split (-> per-hit files). Both for exact and local.

## 8a. Finding: local s1==s2 is a FREE-FIELD artifact (2026-06-12 probe)
Measured (`test_samesite_spindiag_claude.cu`, same-site block `B(t)=[D_m^{-1}]_{(n0,t),(n0,t0)}`,
`max|offdiag|/max|diag|`): U=1 -> 3.6e-16 (spin-diagonal => s1==s2 EXACTLY); Gaussian U -> ~1e-3
(w=0.1/0.3/1.0 give 1.1e-3/2.2e-3/1.6e-3, saturating ~1e-3, NOT growing with the gauge width).
- WHY: at U=1 each site is isotropic (balanced links) so the same-site return paths cancel -> block is
  `a*1 + d*sigma_3` (diagonal) -> sigma_1 B sigma_1 = sigma_2 B sigma_2 -> s1==s2 (signal AND noise).
  A gauge field breaks the per-site isotropy -> small off-diagonal (spin-flip) -> s1 != s2 per config.
- CONSEQUENCE: NO exact saving for interacting (keep s1 and s2 both for an exact result). But s1 ~= s2 to
  ~0.1-0.2%, FAR below the stochastic error (>>1%), so the local spatial current may be approximated
  `G_s^loc = (s1+s2)/2 ~= s1`, dropping the s2 channel solves (~n_sites/pattern, ~1/3 of LOCAL source
  solves, ~9% of total source solves; ~constant in L).  OPTIONAL optimization to revisit if the local-Gs
  noise budget is tight; the dominant cost (exact-K + axial) is untouched.  Default: KEEP all three
  channels (the free-field validation needs s1,s2,s3 separately anyway).

---

## 9. Status (2026-06-12: IMPLEMENTED + free-field validated at 1 hit)
- **Chunks 1-3L DONE** in `jj_corr_dilute_claude.cu` (builds clean; free L=1 run status 0):
  Ch1 header+`seed_from_string`+`corr_dil` tag; Ch1b ylm stripped (`N_ELL=0`); Ch2 hoisted single
  accumulators (Ctp/Csp/Atp/Asp/Jtp/Jsp + local Cs/Js); Ch3 summed-origin SUPERPOSED source + 4-pattern
  spin x e/o DILUTION loop + per-hit string reseed; Ch3L local s1,s2,s3 conn + local disc
  (w_site=dual_areas); Ch4 per-hit write (`corr.<k>.h<h>.h5`, t0s + rng_seed).
- **VALIDATED (free, 1 hit):** exact tp/sp track `corr_deter_exactsum_L1`; superposition wrap confirmed
  (`C[Nt/2-dt] ~ exactsum[dt]`); local s3<0, s1=s2>0 track `corr_deter_local_L1`; s1==s2 to 1e-16 (free
  spin-diagonal block, Sec. 8a).  Parity blocks left UNCONVERTED to superposition (dead in free; PARITY-TODO).
- **Stage 4 DONE:** `run_jj_dilute_claude.sh [GPU NHITS NT0 ENSDIR]` (free, or interacting via --ens-dir).
- **Stage 5 DONE:** `jj_dilute_validate_claude.ipynb` (exact vs exactsum, local vs corr_deter_local, ratios
  ->-1, s1==s2 check, diluted-vs-volume variance at matched cost).  Needs a MULTI-HIT free run for
  jackknife: `bash run_jj_dilute_claude.sh 0 16 2`.
- Side probe: `test_samesite_spindiag_claude.cu` + `tmp_spindiag_claude.sh` (Sec. 8a finding).
- TODOs: multi-hit statistical validation; parity superposition; interacting/mass/R; optional s2~=s1
  saving (Sec. 8a); local ylm + local axial (`jj_stoch_local_impl_plan_claude.md`).
