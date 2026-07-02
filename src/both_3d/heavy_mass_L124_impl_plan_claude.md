# Heavy-mass L=1,2,4 study (Nf=2) — impl plan

**Date:** 2026-06-26. **Goal:** a NEW, separate ensemble study — heavy masses, Nf=2, at L=1,2,4,
gsq=8.0, MPS-packed 3-per-GPU, kmax=320, funded by the **affine** allocation (freed now that the
L2 campaign is complete). Fully isolated from the running diagonal-mass campaign.

## Masses (decided)
Anchor **effective(L=1) = 0.4, 0.8, 1.2**; the operator's L=1 measure factor is
mean_dual_area/mean_ell = 0.9458508 (= 1/1.0572491470487), so the **physical m (R=1, same at every
L)** are effective_L1 x 1.0572491470487:
| physical m (CLI mass_re) | eff L1 | eff L2 (x0.506305) | eff L4 (x0.259021) |
|---|---|---|---|
| 0.4228996588195 | 0.400 | 0.2141 | 0.1095 |
| 0.845799317639  | 0.800 | 0.4282 | 0.2191 |
| 1.2686989764584 | 1.200 | 0.6423 | 0.3286 |
All heavy at every L -> well-conditioned overlap, no Wilson zero-crossing risk (wide k=0.001 not
needed; per-L current params suffice).

## Decisions (NM, 2026-06-26)
1. **Per-L current params**: L1 n=21/k=0.01/nsteps(Nf2)=8 ; L2 n=17/k=0.01/nsteps=6 ;
   L4 n=21/k=0.001/nsteps=8. (All inherit the no-adaptive shared header `overlap_wmass_claude.h`.)
2. **All three L on the affine account** (`--account=affine.lq2_gpu`).
3. **Option (a): SEPARATE heavy binaries** with kmax=320, so the running campaign binaries are
   untouched (esp. L4, whose live light-mass campaign stays capped at 120).

## Binaries
- **L1**: NEW driver `hmc_fermilab_wmass_heavy_claude.cu` = copy of `hmc_fermilab_wmass_claude.cu`
  with **kmax 1200 -> 320** (`:209`). Already n=21/k=0.01/nsteps(Nf2)=8. (Optional: add the
  `# mass_coeff` startup print the L2/L4 drivers have, to confirm effective=0.4/0.8/1.2.)
- **L4**: NEW driver `hmc_fermilab_wmass_L4_heavy_claude.cu` = copy of
  `hmc_fermilab_wmass_L4_claude.cu` with **kmax 120 -> 320** (`:221`). Keeps n=21/k=0.001/
  nsteps(Nf2)=8, `#define GRAD_L4`, N_REFINE=4.
- **L2**: **REUSE the existing** `hmc_fermilab_wmass_L2_claude.o` (already kmax=320,
  n=17/k=0.01/nsteps=6; new heavy masses land in NEW dirs -> no collision, no rebuild needed).
- Build (NM, GPU node): `make -f Makefile.fnal hmc_fermilab_wmass_heavy_claude.o` and
  `... hmc_fermilab_wmass_L4_heavy_claude.o` (generic `%.o:%.cu` rule). Reference the .o from the
  src dir `/project/qed3/qed3/src/both_3d/`.

## Run scripts (live in a dedicated `/lustre2/affine/heavy/` so dirs/logs stay separate)
- **`run_heavy_mps3_claude.sh`** — 3-client MPS run script (generalizes the 2-client
  `run_nf_fermilab_mps_claude.sh`): one MPS daemon, **3 backgrounded clients** (masses M1,M2,M3),
  each with its own `MAX_SEC` 6th arg + `BACKSTOP` timeout, each writing its own checkpoint dir
  (dir3 encodes mRe -> distinct). Header: `--account=affine.lq2_gpu`, `--partition=lq2_gpu`,
  `--gpus=a100:1`, NO `--time` in header (wrapper passes `--time` per L). Reads via `--export`:
  `BIN, M1, M2, M3, Nf, MAX_SEC`. CLI per client: `$BIN 8.0 $Nf 1.0 $Mi 0.0 $MAX_SEC`.
- **`run_wrapper_heavy_claude.sh`** — for L in 1 2 4: pick `BIN` (L1/L4 heavy .o, L2 existing .o)
  and `WALL`/`MAX_SEC` (L1 4h, L2 4h, L4 16h), submit `run_heavy_mps3_claude.sh` with
  `--time=$WALL --job-name=hmc_heavy_L${L}` and the 3 masses, chained `afterany` NCHAIN-deep
  (same `--parsable` PREV[] bookkeeping as the campaign wrapper). Optional `L_LIST` env override.

## Cost / sanity
- L1 tiny (12 sites) -> 320x3 finishes in ~one short job. L2 small. L4 heavy-Nf2 is the cheapest
  L4 corner. 3-client MPS fits easily on one A100 (L4 ~1.3GB/client -> ~4GB of 80GB).
- Startup checks (per client log): `# alat`, admissibility assert PASS, `# delta` (small, heavy
  mass), `# physical_m` = the value above; for L1/L4 confirm effective via the mass_coeff/measure.

## Ordered chunks
1. **Drivers** — create the 2 heavy `.cu` (L1 + L4, kmax=320). Files:
   `hmc_fermilab_wmass_heavy_claude.cu`, `hmc_fermilab_wmass_L4_heavy_claude.cu`.
2. **Scripts** — `run_heavy_mps3_claude.sh` + `run_wrapper_heavy_claude.sh` in `/lustre2/affine/heavy/`.
   Files: those two.
3. **Build + launch handoff** — `/lustre2/affine/heavy/tmp_claude.sh`: build the 2 heavy .o
   (early-exit on fail), then NM launches the wrapper. Dry-run default.

## Open questions
- Walltimes L1/L2/L4 = 4h/4h/16h OK? (L1 could be <1h; 4h is a safe ceiling, graceful-stop handles it.)
- NCHAIN per launch (chain depth)? Default 4 like the campaign.
- `/lustre2/affine/heavy/` subdir OK as the home for this study's dirs+scripts+logs?
