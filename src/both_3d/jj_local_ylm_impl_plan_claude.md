# Impl plan — Ylm-projected tower Eq.(4.35) in jj_local_deter (local op, cont prop, deterministic trace)

Add the spherical-harmonic descendant tower $G^t_{\ell\ell}(t)$ (qed3int_v2-11 Eqs. 4.33-4.35) to
`jj_local_deter_claude.cu` (local $\sigma^a$ current, deterministic trace, continuum propagator P) and
the comparison cells to `checkCFT_claude.ipynb`. Reference for the validated continuum tower:
`jj_cft_ylm_check_claude.{cc,md}` (same propagator, exact Gauss-Legendre integral).

## Physics / goal

The Ylm tower is the spherical-harmonic decomposition of the **temporal** current correlator
$G_t(t;\hat n_1,\hat n_2)$ (= the $\sigma_3$ / "s3" channel between arbitrary site pairs). The paper's
diagonal projection (4.33), using isotropy + Funk-Hecke, is the single Legendre coefficient

$$
g_\ell(t) \equiv G^t_{\ell\ell}(t) = \tfrac12\int_{-1}^{1} G_t(t;x)\,P_\ell(x)\,dx,
\qquad x=\hat n_1\!\cdot\!\hat n_2 .
$$

Prediction (4.35), $\Delta=2$, $C_j=-1/(8\pi^2)$ (from `jj_cft_ylm_check`):

$$
G^t_{11}(t)\to \tfrac{1}{24\pi^2}e^{-2t},\qquad
G^t_{22}(t)\to \tfrac{1}{10\pi^2}e^{-3t},\qquad
G^t_{00}(t)\to 0 .
$$

Convention-free checks: decay rates $(2,3)$ for $\ell=(1,2)$; ratio
$[G^t_{22}e^{3t}]/[G^t_{11}e^{2t}]\to 12/5=2.4$; $G^t_{00}/G^t_{11}\to0$.

## Deterministic lattice estimator (diagonal-m, matches 4.35)

Lattice quadrature of $g_\ell$ ($\int d^2n\to\sum_n A_n$, $A_n=$ `dual_areas[n]`):

$$
g_\ell(t)=\frac{1}{(4\pi)^2}\sum_{n,n'}A_n A_{n'}\,P_\ell(\hat n\!\cdot\!\hat n')\,
          \mathrm{tr}[\sigma_3(n,t_0)\,P\,\sigma_3(n',t)\,P].
$$

Using the addition theorem $P_\ell(\hat n\!\cdot\!\hat n')=\frac{4\pi}{2\ell+1}\sum_{m} Y_{\ell m}(\hat n)Y_{\ell m}(\hat n')$
(real $Y_{\ell m}$, `Ylm_real` in `valence_claude.h`), this FACTORIZES per $m$:

$$
g_\ell(t)=\frac{1}{4\pi}\,\frac{1}{2\ell+1}\sum_{m=-\ell}^{\ell}
          \mathrm{tr}[\Sigma_{\ell m}(t_0)\,P\,\Sigma_{\ell m}(t)\,P],\qquad
\Sigma_{\ell m}(t)=\sum_n A_n\,Y_{\ell m}(\hat n)\,\sigma_3(n,t).
$$

$\Sigma_{\ell m}(t)$ is a sparse matrix (the bare $\sigma_3$ site blocks of `local_W_sigma`, scaled by
$A_n Y_{\ell m}(\hat n)$); the trace reuses `tr_WPWP` unchanged. The accumulator stored is
$C_\ell[dt]=\frac{1}{2\ell+1}\sum_m \mathrm{tr}[\Sigma_{\ell m}(t_0)P\Sigma_{\ell m}(t)P]$; `write_corr`
folds the remaining $1/(4\pi)$ (same as s1/s2/s3). This reproduces `jj_cft_ylm_check`'s $g_\ell$ up to
the 12-site (L=1) quadrature error.

**Normalization choice.** This uses the diagonal-$m$ / Legendre-coefficient $g_\ell$ normalization
(the $1/(2\ell+1)$ + addition theorem), so the output matches Eq.(4.35) directly (ratio 2.4). This is
DIFFERENT from the stochastic `jj_corr_block_t` ylm tower, which uses the outer-product weight
$W_\ell(n)=A_n\sum_m Y_{\ell m}(n)/\kappa$ (no $1/(2\ell+1)$, no diagonal-$m$); that convention gives
$(2\ell+1)g_\ell$ (ratio $5/3\times2.4=4.0$) and agrees with $g_\ell$ only in the continuum (isotropic
kernel). Open question O1 below.

## Files to modify

- `jj_local_deter_claude.cu` — add the ylm tower loop + output (chunk 1).
- `checkCFT_claude.ipynb` — add ylm load + Eq.(4.35) comparison cells (chunk 2).

## Chunk 1 — C++ (jj_local_deter_claude.cu)

Files: `jj_local_deter_claude.cu`

1. After the s1/s2/s3 channel loop, add a `// ---- Ylm tower (Eq. 4.35), temporal sigma_3 channel ----`
   block, `constexpr int L_MAX_YLM=2` (l=0,1,2).
2. Precompute `A_n = base.dual_areas[n]` (already `w_site[n]`) and, per (l,m,n),
   `ylw[l][m][n] = w_site[n]*Ylm_real(l,m,base.sites[n])`.
3. For each `l`, init `Cyl[b][dt]=0`; for each `m in [-l,l]`:
   - build `Sig0[b]` = $\Sigma_{\ell m}(t0_b)$ and `Sigt` = $\Sigma_{\ell m}(t)$ as `Ent` lists:
     per site n push `{off,off,+ylw}` and `{off+1,off+1,-ylw}` (= $\sigma_3$ scaled), `off=Nx*t+NS*n`.
   - `Cyl[b][dt] += tr_WPWP(Sig0[b], Sigt, P)`.
   - after the m-loop divide by `(2l+1)`.
4. Output: `write_corr(h5, "h0/t0_"+b+"/ylm/l"+l+"/Vpp", Cyl[b], false)` and `/Vmm` (conj) if `!parity`.
   Add `ls` already written ({0,1,2}); reuse it.
5. A small file-scope helper `build_Sigma_ylm(en, w_site, sites, l, m, t)` (named static fn, no lambda)
   keeps the loop clean.

Cost (L=1): l=0..2 (9 m), Nt=128, n_t0=2, tr_WPWP O((2 n_sites)^2)=O(576) -> ~1.3M flops. Trivial.

## Chunk 2 — notebook (checkCFT_claude.ipynb)

Files: `checkCFT_claude.ipynb`

1. `load_Vpp(corrdir, "ylm/l"+str(l))` already works (generic key).
2. New cell: load `g0,g1,g2`; plot `g1` vs `e^{-2t}` and `g2` vs `e^{-3t}` (log scale, amplitude matched
   at one point like the s3/sp cells), title "Eq.(4.35) tower".
3. New cell (convention-free): effective decay rates of g1,g2 (should -> 2,3); ratio
   `(g2*exp(3t))/(g1*exp(2t))` -> 2.4; `g0/g1` -> 0. Markdown noting the diagonal-m normalization.

## L=1,2,4 extension (added)

The binary is L-parametrized (`-DN_REFINE_CLI`); no per-L code change. Two additions:
- `load_mat` made memory-lean (scope the real/imag buffers separately) so peak RAM = M + one N^2
  double buffer (~41 GB at L=4) instead of M+real+imag (~55 GB) -> the 26 GB `cont_prop_L4` fits on
  the 62 GB host. L=4 is RAM-tight; run when idle or on a larger node.
- `tmp_claude.sh` loops L in {1,2,4}: compile `jj_local_deter_L<L>.o`, run `--prop-file
  cont_prop_L<L>/Dinv.0.h5 --out-tag cont`. (User must `rm` the stale L=1,L=2 `corr.0.h5` first --
  `complete` flag => skip-on-rerun.)
- Notebook: existing L=1 plots relabeled "cont prop at L=1"; new `## L-convergence` section overlays
  G_t, G_s (same-site -> L-independent, collapse) and the Ylm tower + ratio (off-site quadrature ->
  converges to rates (2,3) and 2.4 as L grows). Missing-L dirs are skipped via try/except.

## Open questions

- **O1 (normalization). RESOLVED = (a) diagonal-$m$ $g_\ell$** (1/(2l+1) sum_m <j_m j_m>; matches
  Eq.4.35, ratio 2.4, reproduces `jj_cft_ylm_check`). User confirmed the structural distinction:
  (a) = sum_m <j_m j_m>, (b) = sum_{m,m'} <j_m j_m'> (off-diag vanish for isotropic G_t).
- **O2 (output key).** `h0/t0_b/ylm/l{l}/{Vpp,Vmm}` (chosen) vs flat `yl{l}` like s{a}. -> plan uses
  `ylm/l{l}`.
- **O3 (channel).** Tower built on the temporal $\sigma_3$ (= G_t) only, per Eqs.(4.33-4.35). Confirm
  (no spatial-current ylm requested).
