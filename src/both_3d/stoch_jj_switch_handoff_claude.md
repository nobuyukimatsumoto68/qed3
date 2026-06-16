# Stochastic JJ — switch-from-deterministic handoff (2026-06-11)

Handoff for the next session: we have finished the **deterministic** free-field validation of the
conserved-current correlators (loc / disp / exact, Eqs. 4.28/4.29/4.36 of **qed3int_v2-13**, source
`qed3int_v2(6)/main.tex`).  The next job is to move to the **stochastic** production code and carry over
the conventions/caveats we settled.  This doc maps the stochastic status and the carry-over caveats.

PDF: `qed3int_v2-13.pdf` ; LaTeX source `qed3int_v2(6)/main.tex`.  Memory: `project_state.md` (top entries),
`project_mrhs_c6.md`, `project_inner_pole_multishift.md`.

---

## 1. Deterministic <-> stochastic correspondence (the key picture)

The **stochastic** program `jj_corr_block_t_claude.cu` computes the **exactly-conserved current K**
correlator stochastically.  It is the random-source, production (interacting + mass + axial) version of
the deterministic **exact** path we just built (`jj_exact_diag_deter_free_claude.cu --sum`,
`corr_deter_exactsum_L<L>`).  They are the **same observable**, Eq.(4.29)/(4.32)/(4.36):

- stoch sp estimator  `(1/4pi) sum_links (link_volume/kappa^2) <K^{nn'}(0) K^{nn'}(t)>`
  = det exactsum sp `(1/4pi) sum_links link_volume <(K/kappa)(0)(K/kappa)(t)>`   (kappa in weight vs in operator).
- stoch tp           `(1/4pi) sum_sites (dual_areas/kappa_t^2) <K^t K^t>`  = det exactsum tp.
- stoch ylm          `W_ell = (A_n/kappa_t) sum_m Y_lm`  = det exactsum Sigma_{l,m} = sum_n A_n Y_lm (K^t/kappa_t).

So the **deterministic `corr_deter_exactsum_L{1,2}` is the ground-truth** the free stochastic run must
reproduce (within stochastic noise).  There is NO stochastic "loc" or "disp" — stochastic = K only.
(loc = on-site sigma, disp = displaced link, exact = full K were three *deterministic* discretizations;
only K is done stochastically.)

---

## 2. CAVEATS to carry over (read before coding the stochastic analysis)

**(A) v2-13 (D-1): the stochastic ratio G_s/G_t = -1 RAW — do NOT divide by (D-1).**
v2-13 defines the spatial correlator as the per-direction average `G^s = (delta^ab - e3 e3)/(D-1) f^{ab}`,
so `G_s/G_t = -1` (sign opposite, equal magnitude).  The `1/(D-1)` adaptation we applied in the
deterministic notebooks is **ONLY for loc**, because loc sums BOTH spatial Pauli channels `s1+s2 = (D-1)g`
(raw ratio -2).  The **link / site point-split currents (disp, exact, and the stochastic K)** are a single
direction per link; the diagonal `<j_e j_e>` on one link is already the per-direction value g
(`e` unit, `<j^a j^b> = g delta^ab`), and the area-weighted link sum averages to g — it does **not** pick
up `(D-1)`.  Verified deterministically: disp/exact **summed** raw `G_s/G_t = -1.07 (L1) / -1.06 (L2)` -> -1
(NOT -2).  ==> The stochastic `w_sp` link sum already gives ratio **-1**; applying `/(D-1)` to it would be
the wrong "-0.5" bug we hit and fixed in the summed notebook.

**(B) Normalization is in the WEIGHT (stoch) or the OPERATOR (det) — equivalent.**
Stoch uses RAW K (kappa-in) with weight `A/kappa^2`; det exactsum uses `K_ov_kappa = K/kappa` with weight
`A`.  `K^2/kappa^2 = (K/kappa)^2`, so identical.  The `_ov_kappa` machinery
(`ConservedCurrent::build_W_ov_kappa`, `insertion_kappa`) is a det convenience; the stoch's `w/kappa^2` is
the same normalization.  Do not double-remove kappa.

**(C) Sign + tower (v2-13).**  `tp` (G_t) Re<0 ; `sp` (G_s) Re>0 ; `G_s/G_t -> -1`.
Ylm tower (Eq.4.36, temporal-current descendants): `g_1 -> -1/3 C_j e^{-2t}`, `g_2 -> -4/5 C_j e^{-3t}`,
`g_0 -> 0` (charge conservation; icosahedral l=6 aliasing, shrinks ~1/L^4).  Convention-free checks:
rates (2,3); `g_2 e^{3t}/g_1 e^{2t} -> 12/5 = 2.4`; amplitude vs same-source G_t: `g_1/G_t -> 1/3`
(since `G_t = -C_j e^{-2t}(1-e^{-t})^{-4}`).

**(D) i-convention is clean in stoch.**  `ConservedCurrent::build_W` multiplies `d_coo_format (= i*C)` by
`i/lambda_M = -C/lambda_M` (the i is removed); the block kernel `build_WH` mirrors it; weights are pure
real.  No floating i (the audit found only the old hand-rolled disp had one, since fixed).

**(E) CSR-build perf fix is live.**  `jj_corr_block_t_claude.cu` (and all `_claude` consumers) now
`#include "sparse_dirac_claude.h"` — O(len) bucketing of the CSR index build (was O(N*len) at startup);
byte-identical (verified, `-DCSR_VERIFY`).  No action; just don't revert.

---

## 3. Stochastic code map

**Production program:** `jj_corr_block_t_claude.cu`  (C6f-c t-blocked-sink variant; the one to use).
Copies: `jj_corr_block_t_fermilab_claude.cu` (abs geometry path, excluded from make), `_L2`/`_L4` per-L
copies, `jj_corr_L4_claude.cu` (non-blocked L4 fallback, fits memory but slow).

**K kernel:**
- `ConservedCurrent<Fermion,Gauge> kop(Dm)` (`conserved_current_claude.h`) — K is mass-independent;
  multishift `apply_k_ms` via `operator()` (inherited through `op_K.from_cpu`).
- `ConservedCurrentBlockT` (`conserved_current_block_claude.h`) — the SINK lever: `apply_k_block_t` /
  `apply_k_dag_block_t` compute `K(t, fixed) xi` for ALL t in one pass (uses `kop.build_W` + own
  `build_WH`).  C6f t-blocking gave the ~3x jj speedup (see `project_mrhs_c6.md`).

**What one run computes per config (UNIFIED, Sec. 3.8):** disc + connected, vector + axial, tp / sp / ylm,
sharing `phi' = D_m^{-1} eta` and the `K(a,t) phi'` sink applies.
- source leg `psi_a(t0) = D^{-dag} K^dag(a,t0) eta` ; sink `kphi = K(a,t) phi'`.
- disc traces `J(t) = sum_a w_a eta^dag K(a,t) phi'` ride the sink pass (no extra K applies); RAW (no 1/4pi).
- connected `Vpp` (and `Vmm`) folded by `1/(4pi)` in `write_corr`.

**Weights (geometry; `jj_corr_block_t_claude.cu:332-351`):**
- `w_tp[n]  = dual_areas[n] / kappa_t(n)^2`            (Eq. 4.32, temporal)
- `w_sp[il] = link_volume[il] / bd.kappa[il]^2`        (Eq. 4.29, spatial)
- `W_ell[l][n] = (dual_areas[n] / kappa_t(n)) * sum_{m=-l}^{l} Y_lm(n^)`  (Eq. 4.36, ylm, kappa^1)

**Mass cases (valence via `--mass-re/--mass-im`):**
- vector `(++)` both legs use `D_m = D_ov + m` (Eq. 3.60).
- parity `m = i m_P` (purely imaginary): dagger leg uses `tilde D_{m_P} = D_ov + m_P/(1-m_P)` (Eq. 3.63);
  needs its own `tilde D^{-1} eta` solve for the parity disc.
- flavor `m = m_F` (purely real): axial uses massless `D_ov` legs (no conserved axial current).
- axial: connected `C_{A+-}` only (`Apm`), GW factor `op_oneMinusD = (1 - D_ov)`.

**Modes / CLI:** `--gsq --Nf --nu0 --nu1 --mass-re --mass-im --ens-dir(OMIT=free U=1) --nhits --n-t0
--ninter --kmin --kmax`.  free-field = omit `--ens-dir` (the CFT check).

**Output:** `data_<ESNID>/corr_nt0<N>_nhits<H>/corr.<k>.h<h>.h5` — **ONE FILE PER HIT** (atomic .tmp+rename,
own `complete` sentinel, resume-skips finished hits); data under `h0/`.  Keys:
`h0/t0_<b>/{tp,sp}/{Vpp,Vmm}`, `h0/t0_<b>/ylm/l<l>/{Vpp,Vmm}`, axial `.../Apm`, disc `h0/disc/...`.
`Vmm = conj(Vpp)` for massless/m_F (`write_corr_conj`).  ESNID: free=`free_vmRe..vmIm..`, else ensemble.

**Run scripts:** `run_jj_analysis_claude.sh` (--gpu, --nhits, --n-t0, tmax), `run_obs_local_claude.sh` /
`run_obs_fermilab_claude.sh` (sea ensembles, per-GPU Nf split), `run_exactjj_free_claude.sh` (free L
copies).

**Analysis notebooks:** `jj_corr_massless_claude.ipynb`, `jj_corr_mF_claude.ipynb`, `jj_corr_mP_claude.ipynb`
(load-only; assert mass pattern; conn tp/sp/ylm vector Vpp/Vmm + axial Apm, disc, sum `-conn+disc`
Eq.3.39; `resample` knob: free=jackknife HITS, interacting=CONFIGS).  R folded into the mP notebook.

---

## 4. Validation plan (do this first)

1. Run the stochastic in **FREE mode** (omit `--ens-dir`) at L=1,2 with enough hits, and compare its
   `tp/sp/ylm` + ratio against the deterministic **`corr_deter_exactsum_L{1,2}`** (ground truth).  They are
   the same observable -> must agree within stochastic error.  Expect `G_s/G_t -> -1`, `tp<0/sp>0`, ylm
   rates (2,3), ratio 2.4, `g1/Gt -> 1/3`, `g0 -> 0`.  **No `/(D-1)` on the stochastic** (caveat A).
2. Confirm the connected-vs-disc sign convention matches (`physical = -C_conn + C_disc`; det exactsum
   writes the connected `Vpp` directly; free disc ~ 0 for gauge-fixed U=1).
3. Then move to interacting (sea `--ens-dir`), mass (m_F / m_P), axial, and reweighting R.

---

## 5. Known open issues / TODO (from memory + this session)

- **Free disc/sum errors need raw-J hit-resampling** — `jj_disc_postproc` hit-collapses the disc trace, so
  free disc/sum jackknife over hits is not yet wired (conn tp is fine: 64 hit-samples ~6-7 sigma).
- **Jackknife resample axis**: free => HITS (gauge fixed, 1 config), interacting => CONFIGS.
- **Reweighting R** (parity sea): `reweighting_R_claude.cu` -> `R/R.<k>.h5`; folded into `jj_corr_mP` nb.
- **Deterministic exactsum L=2** is ~1.6 day (build-use-discard, ~1.74M op_K solves); L=1 ~2.8 h.  Useful
  as the L=2 cross-check if patience allows; otherwise validate at L=1 and trust the L=2 stochastic.
- **exact `--sum` ylm**: the det exactsum now also emits the ylm tower (folded into `--sum`); restart the
  exact run on the ylm-augmented binary to get it (see `project_state.md` top entry).

---

## 6. One-line orientation for the next session

> The stochastic `jj_corr_block_t` already computes Eq.(4.29/4.32/4.36) for the exact current K with the
> RAW-current + `A/kappa^2`-weight form; it equals the deterministic `corr_deter_exactsum`.  Validate free
> stochastic vs that ground truth at L=1,2, **remembering G_s/G_t = -1 with NO `/(D-1)` for the link
> current** (the `/(D-1)` was a loc-only thing), then go interacting + mass + axial + R.
