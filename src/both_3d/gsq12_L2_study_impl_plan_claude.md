# gsq=12, L=2 study -- implementation plan

Created 2026-07-02. A new coupling point for the diagonal measure-weighted L=2 ensembles
(sibling of the gsq=8.0 campaign [[l2l4-refinement-run]] and the heavy study
[[heavy-mass-l124-study]]).

## Goal / physics

Generate L=2 (radial-quantization S2xR, $N_\text{refine}=2$) HMC ensembles at **stronger
coupling $g^2=12$** (campaign is $g^2=8$), diagonal measure-weighted overlap mass
$m_L = m\,A_y/\bar a_s$ (arXiv:2510.03085). Same 7 physical masses as the campaign+heavy
sets, all at $N_f=2,4,6$, cap 320 configs each.

**Zolotarev $n=21$, window $k=0.001$** (NM decision). Rationale: at higher $g^2$ the gauge
fields are rougher, so the mass-independent Wilson kernel $D_W^\dagger D_W$ can have lower
low-modes -> a wider fixed window + more poles guards the sign approximation against the
zero-crossing spikes seen at L=4 (same fix as [[l2l4-refinement-run]]; Zolotarev rational
approx of $\text{sign}(H_W)$, A.D.Kennedy hep-lat/0402038). The campaign L2 used $n=17/k=0.01$
(fine at $g^2=8$); this study overrides to the safe wide window.

## Masses (physical $m$, R=1 units -- SAME at every L; the factor $A_y/\bar a_s$ is L-dependent)

- **4 light -> qed3 allocation** (2 MPS pairs, 2 clients/GPU, like the campaign):
  pairA = {0.0105724914705, 0.2114498294097}, pairB = {0.0528624573524, 0.1057249147049}
- **3 heavy -> affine allocation** (1 MPS triple, 3 clients/GPU, like the heavy study):
  {0.4228996588195, 0.845799317639, 1.2686989764584}  (effective_L1 = 0.4/0.8/1.2)

NB account split (NM, final): **light -> qed3, heavy -> affine**. $N_f=\{2,4,6\}$ for all 7.
=> 4 light x 3 Nf = 12 light streams (6 pair-jobs), 3 heavy x 3 Nf = 9 heavy streams (3 triple-jobs).

## Params

gsq=12.0 (CLI arg), nu0=1.0, Nt=128, at=0.2, L=N_refine=2 (all unchanged from campaign L2 driver).
**kmax=320** (already the L2 driver value). **nsteps=6** all Nf (current L2 driver value; watch
acceptance -- higher gsq may need more steps; bump only if acc drops, do NOT preempt). tmax=1.0.

## Files

### 1. New binary source (one binary serves BOTH light & heavy; gsq is a CLI arg)
`~/project_qed3/qed3/src/both_3d/hmc_fermilab_wmass_L2_gsq12_claude.cu`
= copy of `hmc_fermilab_wmass_L2_claude.cu`, ONE change at the Zolotarev line (~:187):
  keep `// Fermion D(DW, mass, 17);` commented, add `Fermion D(DW, mass, 21, 0.001);`.
Everything else identical (kmax=320, nsteps=6, NPARALLEL_DUPDATE=1 already for MPS).
Build: `make -f Makefile.fnal hmc_fermilab_wmass_L2_gsq12_claude.o`.

### 2. Run dirs (checkpoints + logs land in cwd = submit dir; dir3 also encodes gsq so no collision)
- `/lustre2/qed3/gsq12/`   (light)
- `/lustre2/affine/gsq12/` (heavy)

### 3. Run scripts
- `/lustre2/qed3/gsq12/run_gsq12_mps_claude.sh` = copy of the campaign 2-client MPS script with
  `--account=qed3.lq2_gpu`, `GSQ=12.0`, BIN=gsq12 .o. --time=4:00:00, WALL_SEC=14400.
- `/lustre2/affine/gsq12/run_gsq12_mps3_claude.sh` = copy of the heavy 3-client MPS script with
  `--account=affine.lq2_gpu`, `GSQ=12.0`, BIN=gsq12 .o. --time=4:00:00, WALL_SEC=14400.

### 4. Wrappers (afterany-chained, --parsable PREV[] scheme; jobnames encode gsq12)
- `/lustre2/qed3/gsq12/run_wrapper_gsq12_light_claude.sh` -- loops Nf={2,4,6} x tag={A,B},
  2-client pairs, jobname `hmc_gsq12_Nf{n}_L2_pair{A,B}`.
- `/lustre2/affine/gsq12/run_wrapper_gsq12_heavy_claude.sh` -- loops Nf={2,4,6}, 3-client triple,
  jobname `hmc_gsq12_Nf{n}_L2_heavy`.

### 5. Handoff `/lustre2/affine/gsq12/tmp_claude.sh` (DRY-RUN default, --apply [NCHAIN])
Build the shared .o once (early-exit on fail), then run BOTH wrappers (light submits to qed3,
heavy to affine -- account is baked into each run script, so one node submits both). Tees a log.

## Chunks

1. **Binary**: create+build `hmc_fermilab_wmass_L2_gsq12_claude.cu` (n=21/k=0.001). Verify it
   compiles; NM confirms `# delta` + `# gsq = 12` at a quick startup check.
2. **Light side (qed3)**: dir + `run_gsq12_mps_claude.sh` + `run_wrapper_gsq12_light_claude.sh`.
3. **Heavy side (affine)**: dir + `run_gsq12_mps3_claude.sh` + `run_wrapper_gsq12_heavy_claude.sh`.
4. **Handoff + launch**: `tmp_claude.sh`; NM runs `--apply`.
5. **Skill + memory**: extend `check-ensembles` to audit the gsq12 study; add a memory file.

## Open questions (resolved)

- Accounts: light->qed3, heavy->affine (NM, final). Masses: 4 campaign light + 3 heavy (all 7).
  Nf=2,4,6. Zolotarev n=21/k=0.001. nsteps: keep 6 (watch acc).
