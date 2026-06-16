# Per-mass validation notebooks -- implementation plan

## Goal
Reorganize the three observable-split validation notebooks (`jj_dilute_validate_claude.ipynb` vector,
`jj_dilute_validate_axial_claude.ipynb` axial, `condensate_validate_claude.ipynb`) into **three
mass-split notebooks**, each self-contained for one mass and containing **all three observable types**
(vector, axial, condensate):

- `jj_validate_m0_claude.ipynb`  -- massless ($m=0$)
- `jj_validate_mF_claude.ipynb`  -- flavor-breaking real mass $m_F=0.1$
- `jj_validate_mP_claude.ipynb`  -- parity-breaking imaginary mass $m_P=0.1i$

The three are structurally **identical** (maximally parallel); they differ only by the `MASS` value in
the config cell and by cells that auto-gate on `MASS` (see Differences).

## Uniform conventions (user directives -- apply in ALL notebooks)
1. **Folding: $N_t/2$-shift translation average only (`half_shift_avg`); NO four-fold averaging anywhere.**
   Every correlator (vector exact tp/sp, vector local s1/s2/s3, axial) is processed identically:
   `half_shift_avg(C)` then plotted on the independent region $0<t<N_t/2$.  The single-origin determ refs
   are put in master-field form via `mirror()` then `half_shift_avg` (a no-op once $N_t/2$-periodic).
2. **Re and Im in separate cells** for every correlator (one concept per cell).  Im is ~noise (around 0)
   for $m_0$/$m_F$ and physical for $m_P$ -- plotted in all cases as a sanity check.
2b. **Signed log axis**: each Re/Im cell is a LOG plot showing BOTH sign branches in the same axes -- the
   positive points as filled markers, $|$negative$|$ as open markers; the determ curve splits into a `+`
   (solid) and `-` (dotted) branch.  (User: "draw +- branches in the same plot.")  Ratios/variance/`s1`-vs-`s2`
   keep their own scales (linear ratios; log magnitudes); condensate stays Re/Im point markers.
3. Errors = jackknife over hits (`jk_cplx` for complex mean + re/im errors).
4. Figures -> `figs/<name>_<MLAB>.pdf`; titles carry `, $m=<MTIT>$`.  `MLAB`/`MTIT` keyed by `MASS`.

## Determ ground-truth normalization (per mass, gated)
- vector exact-K: $m_0$ uses `corr_deter_exactsum_L1` ($\times1$); $m_F$/$m_P$ use `corr_deter_exact1_L1`
  ($\times 4\pi$, since exactsum $=4\pi\cdot$ exact1).
- vector local: `corr_deter_local_L1` (direct, summed norm).
- axial exact-K: `corr_deter_exact1_axial_L1` with the matching constant $c=1/4\pi$ (mass-independent).
- axial local: `corr_deter_local_axial_L1` (direct).
- condensate: dense `condensate_deter_L1/cond.h5` (direct, $/V_{st}$).

## Differences (auto-gated on MASS; present in all three for parallel structure)
- **$V_{--}$ (Vmm) channel**: independent only for $m_P$ ($m_0$/$m_F$: $V_{--}=\overline{V_{++}}$ exactly).
  The four Vmm cells (tp/sp x Re/Im) `print` a note and skip for $m_0$/$m_F$.
- **Variance (diluted vs volume `corr_nt02_nhits64`)**: volume run exists only for $m_0$; cell auto-skips
  (try/except) for $m_F$/$m_P$.
- **CFT-shape overlay**: OPEN QUESTION (see below).

## Cell inventory (identical across the three; ~32 cells)
- title (md), imports, config (MASS), helpers (loaders + jk_cplx + half_shift_avg + mirror + DET readers +
  prep_v vector + prep_ax axial + condensate loaders + CFAC)
- VECTOR (md): exact tp V++ Re, tp V++ Im, sp V++ Re, sp V++ Im, exact ratio Gs/Gt;
  exact tp V-- Re, tp V-- Im, sp V-- Re, sp V-- Im (gated mP);
  local s3 Re, s3 Im, (s1+s2)/2 Re, (s1+s2)/2 Im, local ratio, s1-vs-s2; variance (gated m0); vector summary
- AXIAL (md): tp Re, tp Im, sp Re, sp Im, local s3 Re, local s3 Im, local ratio, axial summary
- CONDENSATE (md): dense+stoch table, sigma_PS (Re/Im), sigma_FS (Re/Im)

## Build method (OPEN QUESTION)
Recommended: a single Python builder `build_validate_nbs_claude.py` that holds the cell list once and emits
all three `.ipynb` (substituting `MASS`) -- DRY, guarantees the three stay parallel, easy to regenerate.
Alternative: hand-build one notebook with NotebookEdit, copy to the other two, retarget `MASS`.

## Resolved decisions (user, 2026-06-15)
1. **Build method**: Python builder `build_validate_nbs_claude.py` (DRY single source of truth).
2. **CFT-shape overlay**: DROPPED everywhere -- compare diluted vs the determ ground truth only.
3. **Old notebooks**: left in place; user will move them to `saved_ipynb/` themselves.  I do NOT delete.
4. **Vmm cells**: `mP` notebook ONLY (4 cells tp/sp x Re/Im).  `m0`/`mF` omit them ($V_{--}=\overline{V_{++}}$).

So the only structural difference between the three: the `mP` notebook has 4 extra `Vmm` cells; otherwise
identical (config `MASS` line aside).  Local-current cells plot the `Vpp` channel only (s1/s2/s3); the
two-channel emphasis is on the conserved exact-K current.  Variance cell (diluted vs volume) is included in
all three but auto-skips where the volume run is absent ($m_F$/$m_P$).  All axes linear (sign-safe);
`s1`-vs-`s2` magnitude plot stays log.
