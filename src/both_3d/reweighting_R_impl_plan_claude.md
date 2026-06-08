# Reweighting factor R -- implementation plan

## Goal

Compute the parity reweighting factor R (PDF `qed3int_v2-10.pdf`, Eq. 2.5) per gauge config, for a
parity mass `mP` read from the CLI.  R is stage 3 of the global pipeline
(`jj_analysis_global_procedure_claude.md`); it multiplies the measured correlators in the notebook.

$$
R \;=\; \left[\frac{\det\big((1-m_P)D_\text{ov}+m_P\big)}{\det\big(D_\text{ov}+m_P\big)}\right]^{\!*}
   \;=\; \left(\prod_i \frac{(1-m_P)\lambda_i+m_P}{\lambda_i+m_P}\right)^{\!*},
$$

where $\{\lambda_i\}$ are the eigenvalues of the **massless** overlap $D_\text{ov}$ (PDF: "calculated
in our case with dense-matrix manipulation").  $R\to1$ as $m_P\to0$.

## Method (mimic `eig_wmass_claude.cu`)

`eig_wmass_claude.cu` already does exactly the hard part:
1. Build the dense $N\times N$ matrix of the operator by applying it to each unit vector
   (`Op.from_cpu<N>(Dxi, e_i)` -> column $i$), `N = Comp::N = 3072` at L=1.
2. All eigenvalues via cuSOLVER `cusolverDnXgeev` (NOVECTOR), `CUDA_C_64F`.

For R we use the **massless** overlap (`Overlap<WilsonDirac>`, i.e. the `NO_MASS` path) so the
spectrum is $\{\lambda_i\}$ of $D_\text{ov}$; `mP` enters ONLY in the post-formula above.

**Better way than recomputing a determinant ratio per mP:** the eigenvalues are mP-INDEPENDENT.
So compute/diagonalize ONCE per config (the expensive step: 3072 overlap applies + one $3072^2$ geev),
**store the eigenvalues**, and evaluate R for the CLI mP (and any other mP later) as a trivial product.
This is the intended dense-spectrum route; no stochastic det-ratio estimator needed at L=1.

## Differences from `eig_wmass_claude.cu` (what to add)

- Use massless `Overlap` (define `NO_MASS`); spectrum independent of valence mass.
- **Read gauge configs**: loop `ckpoint_lat.k` in `--ens-dir` (mirror `jj_corr_claude.cu` config loop +
  `complete`-sentinel resume skip), one R per config.  (`eig_wmass` uses a single default U.)
- After geev, compute `R = conj( prod_i ((1-mP)*lam_i + mP)/(lam_i + mP) )` for the CLI mP.
  Accumulate the product in log space (`sum log`) for numerical safety, then `exp`.
- **Output** per config: the eigenvalues `lam` (real,imag, length N) + `R` (real,imag) + `mP` + the
  `complete` sentinel.  Storing eigenvalues lets the notebook reweight for ANY mP without rerunning.

## Files

- NEW `reweighting_R_claude.cu` -- the R/spectrum program (mimics `eig_wmass_claude.cu`).
- NEW `run_reweighting_R_claude.sh` -- driver (loop ensemble configs; CLI mP; GPU pick).
- This plan.

## Ordered chunks

1. **Spectrum core.** Copy the `eig_wmass` dense-build + geev into `reweighting_R_claude.cu`, massless
   Overlap, get `{lam_i}` for one config.  Files: `reweighting_R_claude.cu`.
2. **Config loop + R + output.** ens-dir config loop (+ resume sentinel), the R product for CLI mP,
   write `R.<k>.h5` (eigenvalues + R + mP).  Files: same.
3. **Run script.** `run_reweighting_R_claude.sh` (ens-dir, mP, GPU).  Files: that script.

## Output layout (RESOLVED -- same obs dir as the correlators)

- **ESNID = `ens + '_vmRe<mass_re>vmIm<mass_im>'`**, identical to `jj_corr_claude.cu`; `mP = Complex(mass_re,
  mass_im)` from the CLI (pure parity => `mass_re=0`, `mass_im=mP`).
- Output goes into the SAME obs directory `data_<ESNID>/` as `corr_*`/`disc_proced`, in a subdir
  `data_<ESNID>/R/R.<k>.h5` (config index `k` = the `ckpoint_lat.k` / `corr.<k>.h5` index).
- Per file: `lam/{real,imag}` (len N eigenvalues of $D_\text{ov}$, mP-independent), `R/{real,imag}`,
  `mP/{re,imag}`, `complete` sentinel (written LAST; resume-skip).  Storing `lam` lets the notebook
  reweight for ANY mP without rerunning the diagonalization.
- **Config source:** interacting => loop `--ens-dir`'s `ckpoint_lat.k`; **free** => the single free config
  only (`k=1`, `kmax=0`, U=1; $R\to1$, a useful check).

## Resolved decisions

1. `mP = Complex(mass_re, mass_im)` (CLI, same convention as jj_corr / eig_wmass).
2. R per config; free = single free config; output in `data_<ESNID>/R/`.
3. Store eigenvalues + R per config (reweight any mP in post).

## Open question

4. **Notebook combination:** reweighted observable $\langle O\rangle_R=\langle O\,R\rangle/\langle R\rangle$
   over configs (program just emits R per config) -- confirm this is the intended downstream use.
