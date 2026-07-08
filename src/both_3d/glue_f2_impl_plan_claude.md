# Implementation plan: $F_{\mu\nu}^2$ action-density glueball ($0^{++}$)

## Goal / physics

Measure the gauge action density $\tfrac14\sqrt{g}\,F_{\mu\nu}F^{\mu\nu}$ on the cylinder $S^2\times\mathbb{R}$
as a scalar ($0^{++}$-like) glueball interpolator, projected onto the $Y_{lm}$ tower, and extract the
lowest scalar mass with the same multi-flow + GEVP + vacuum-subtraction pipeline already used for the
linear $F_{12}$ operator in `glue2_msm_claude.cu` / `glue_ylms3_claude.ipynb`.

**Source of the operator definition.** N. Matsumoto et al., `qed3_v2-6.pdf`, Sec. IV, Eqs. (IV.24)-(IV.35).
The lattice Wilson gauge action is the sum over prisms $(\triangle,t)$ of the discretized $F^2$ density.
Per timeslice $s$ the local density is exactly the integrand already coded in `includes/action_ext.h`
(non-compact branch, `is_compact=false`, which is the production setting, `glue2.cu:45`):

$$
\rho(s)\;=\;\underbrace{\sum_{\text{faces }i} \tfrac12\,\beta_s[i]\,\theta_{\rm sp}(s,i)^2}_{\text{magnetic } B^2\ (F_{12}^2)}
\;+\;\underbrace{\sum_{\text{links }il} \tfrac12\,\beta_t[il]\,\theta_{\rm tmp}(s,il)^2}_{\text{electric } E^2\ (F_{0i}^2)}
$$

with `beta_s[i]` (per spatial face) and `beta_t[il]` (per link) the SIMULATED couplings from
`U1WilsonExt::set_beta()` (`action_ext.h:49` spatial, `action_ext.h:65` temporal). $\theta_{\rm sp}$ is the
spatial plaquette angle (`plaquette_angle(s,face)`), $\theta_{\rm tmp}$ the temporal plaquette angle
(`plaquette_angle(s,BaseLink)`). Reusing the action's own `beta_s`/`beta_t` guarantees the measured
density is exactly the local integrand of the simulated action.

The $Y_{lm}$-projected interpolators (mirroring `plaquette_angle_avg_Ylm_real` = face centroid, and
`plaquette_angle_avg_temporal_Ylm_real` = link midpoint, in `gauge_ext.h`):

$$
O^{B^2}_{lm}(s)=\sum_{\text{faces }i} Y_{lm}(r_i)\,\tfrac12\beta_s[i]\,\theta_{\rm sp}(s,i)^2,
\qquad
O^{E^2}_{lm}(s)=\sum_{\text{links }il} Y_{lm}(r_{il})\,\tfrac12\beta_t[il]\,\theta_{\rm tmp}(s,il)^2 .
$$

Only $l=0$ is physically of interest, but the full tower $l=0,1,2,3$ (16 channels, `lm_set` in
`glue2_msm_claude.cu`) is measured, as for $F_{12}$.

## Key differences from the existing $F_{12}$ glue code

1. **Quadratic, scalar, nonzero VEV.** Unlike the linear pseudoscalar $F_{12}$, $F^2$ has a large VEV
   ($\langle\rho\rangle\sim$ mean action density). Vacuum subtraction $C_{ij}=\langle O_iO_j\rangle-\langle O_i\rangle\langle O_j\rangle$
   is essential; the $l{=}0\times l{=}0$ block is dominated by the disconnected constant. The pipeline
   already does this (`vacuum_sub`, notebook cell 27; one-point file `F.k`).
2. **Magnetic + electric.** Two physically distinct pieces ($B^2$, $E^2$) that overlap the same scalar
   state with different weights -> natural variational partners (see Q1).
3. **Signal window is short, $t\approx1$–$4$** (`tcut=5` in the notebook, $a_t=0.2$). GEVP's role is
   small-$t$ excited-state removal so an effective-mass read at $t=1,2$ is trustworthy; no long plateau
   expected. Smearing is spatial-only (`Flow::get_spatial`), orthogonal to the $t$-decay, so it does not
   eat the $t$-window; it only trades ground-state overlap vs. connected-signal magnitude. Multiple flow
   times enter the GEVP basis; the right amount is found by trial and error.

## Parameters (from production, `jj_ensembles_claude.txt`)

- $a_t=0.2$, $Nt=128$, main coupling $g^2=8$ (SEA, mass-free), refinement $L\in\{1,2,4\}$.
- $L$ is compile-time via `-DN_REFINE_CLI` guard (same convention as the L2/L4 stoch drivers).
- Configs read from `Nf{Nf}_gsq{gsq}at0.2nu0{nu0}nt128L{L}/ckpoint_lat.k` (interacting) branch.

## Files to modify / create

- `includes/action_ext_claude.h` (NEW; copy of `action_ext.h` + two added members
  `density_Ylm_spatial(U,s,ell,em)` and `density_Ylm_temporal(U,s,ell,em)`). Non-destructive: original
  `action_ext.h` untouched; only `glue_f2_claude.cu` includes the `_claude` copy.
- `glue_f2_claude.cu` (NEW; copy of `glue2_msm_claude.cu`) — swap the `obs[op][t]` fill from the linear
  `plaquette_angle_avg_Ylm_real` to the $B^2$/$E^2$ action-density channels; add `N_REFINE_CLI` guard;
  output dirs `data_Nf..._f2/` (distinct from `_msm`).
- `glue_f2_claude.ipynb` (NEW) — analysis: opset with the $B^2$/$E^2$ channel layout, vacuum subtraction,
  GEVP, multi-flow, effective mass at $t=1$–$4$.
- `tmp_claude.sh` + `glue_f2_smoke_claude.log` (handoff build+smoke run for the user).

## Ordered chunks

### Chunk 1 — action-density operator methods
Add to `includes/action_ext_claude.h` two members on `U1WilsonExt` that reuse `beta_s`/`beta_t` and the
face-centroid / link-midpoint `Y_{lm}` projection (copied from `gauge_ext.h` patterns), returning
$O^{B^2}_{lm}(s)$ and $O^{E^2}_{lm}(s)$ for one $(s,\ell,m)$.
Files: `includes/action_ext_claude.h`.

### Chunk 2 — driver
Copy `glue2_msm_claude.cu` -> `glue_f2_claude.cu`. Replace the `obs` fill loop (`glue2_msm_claude.cu:298-312`)
with the action-density channels (magnetic and/or electric per Q1). Include `action_ext_claude.h`. Add the
`#ifndef N_REFINE_CLI` guard. Change `dir4` to a `_f2` suffix. Keep `F_corr.k` / `F.k` output format
unchanged so the notebook pipeline is reused.
Files: `glue_f2_claude.cu`.

### Chunk 3 — analysis notebook
New `glue_f2_claude.ipynb` mirroring `glue_ylms3_claude.ipynb` (first creation): opset for the $B^2$/$E^2$
layout, `vacuum_sub`, GEVP with $t_0$ metric, effective mass $-\log\lambda/a_t$, multi-flow basis, folding.
Files: `glue_f2_claude.ipynb`.

### Chunk 4 — build + smoke + flow tuning
Handoff `tmp_claude.sh`: build `glue_f2_claude.cu` (L=1), run on one ensemble's first few configs, verify
`F_corr`/`F` written, VEV-subtracted effective mass finite at $t=1$–$4$. Then trial-and-error over flow
times / N_FLOW.
Files: `tmp_claude.sh`, `glue_f2_smoke_claude.log`.

## Resolved decisions

- **Q1 — SEPARATE $B^2$ and $E^2$ channels** (`NCH=2`, `nops = N_FLOW*NCH*n_lm = 96`, op index
  `iflow*(NCH*n_lm) + ich*n_lm + ilm`, ch 0 = $B^2$, ch 1 = $E^2$). Chosen for the richer variational
  basis; full sum reconstructable in analysis.
- **Q3 — NEW `glue_f2_claude.ipynb`** (first creation; layout differs from the linear code).
- **Q2 — smoke on Nf2 gsq8 L1**; production set (Nf/L) to be confirmed after the smoke run.

## Status (all chunks implemented; build/smoke pending user run)

- Chunk 1 DONE: `includes/action_ext_claude.h` = `action_ext.h` copy + `density_Ylm_spatial` ($B^2$) and
  `density_Ylm_temporal` ($E^2$), reusing `beta_s`/`beta_t`. Original `action_ext.h` untouched.
- Chunk 2 DONE: `glue_f2_claude.cu` (copy of `glue2_msm_claude.cu`): includes `action_ext_claude.h`,
  `N_REFINE_CLI` guard, `NCH=2` action-density `obs` fill, `_f2` output dir, optional 4th CLI arg
  `kmax_run` (cap config count for smoke / flow tuning).
- Chunk 3 DONE: `glue_f2_claude.ipynb` — load `_f2` output, `vacuum_sub`, $l=0$ GEVP over the
  $\{B^2,E^2\}\times$flow sub-basis, single-op effective-mass diagnostic.
- Chunk 4 handoff: `tmp_glue_f2_claude.sh` (build + smoke run, tees to `glue_f2_smoke_claude.log`).
  User runs: `bash tmp_glue_f2_claude.sh` (default L=1 Nf2 gsq8, 5 configs).

## Session 2 update (2026-07-04) — production + HDF5 + C++ GEVP

- Scope grew to BOTH operators on massless L1 gsq8 Nf2/4/6, `kmin=100 stride=1` (all configs).
- `glue2_msm_claude.cu` (linear F_12): `ell=0` DROPPED (monopole flux, not a glueball) -> n_lm=15, nops=45.
- Both drivers: added `kmin/stride/kmax_run` CLI args; REMOVED the k-loop OMP pragma (parallelism is
  ensemble-level, one process per Nf); OUTPUT switched to per-config **HDF5** `F_corr.<k>.h5`
  (datasets `F_corr` Nt x nops^2, `F` nops, `complete` flag; resume via complete-flag + `filesystem::exists`).
- Analysis is C++ (user chose over notebook): `glue_gevp_analysis_claude.cu` — host g++ + Eigen + HighFive,
  streams h5 + per-bin jackknife, faithful glue_ylms3 GEVP (fold/vacuum_sub/zero_noise/t0-metric),
  `inv_sqrt_sym(rtol)` regularization fixes the NaN from the singular 96-op metric. Writes gnuplot `.dat`.
- Handoffs: `run_glue_massless_L1_claude.sh` (measurement, USER runs), `glue_gevp_run_claude.sh` (analysis).
- h5 build flags: `-I/mnt/hdd_barracuda/opt/highfive/include/ -I.../hdf5-2.1.0/include/ -L.../hdf5-2.1.0/lib/ -lhdf5`.
- STATUS: user re-running Step 1 measurement to h5. Then Step 2 analysis -> masses. Perf: text 83GB/ens
  (~17 min I/O) -> h5 (~5 min). Larger-L improvement = next WIP.

## Next (pending)

1. User runs `tmp_glue_f2_claude.sh`; confirm build + `F_corr`/`F` output, VEV-subtracted effective mass
   finite at $t=1$–$4$.
2. Flow-time trial-and-error (`N_FLOW`, `FLOW_INCR` in the driver).
3. Confirm production set (Nf, L) and run full ensembles.
