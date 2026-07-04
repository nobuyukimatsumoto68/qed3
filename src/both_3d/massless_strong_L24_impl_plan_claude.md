# Massless strong-coupling L=2,4 study (gsq=12,16) -- impl plan

Created 2026-07-03. NEW **highest-priority** study (NM): massless (m=0.0) ensembles at strong
coupling, to sit alongside the existing massive campaign / gsq12 / heavy studies. Kill the other
running jobs and focus GPUs here.

## Goal / matrix

Massless overlap ($m=0$) sea ensembles, diagonal measure-weighted scheme (m_L = m*A_y/abar_s = 0):
**Nf{2,4,6} x gsq{12.0,16.0} x L{2,4}** = **12 ensembles**. (L=6 requested but n6 geometry does
not exist -- only n1/2/4/8; NM: do L=2,4 only for now.)

- Caps: **L2 = 620**, **L4 = 320**.
- Zolotarev **n=21, k=0.001** (most conservative from the massive runs; L4/gsq12/heavy all use it).
- nsteps (most conservative seen in the massive runs): **L2 = 8** (all Nf); **L4 = 8/8/10** (Nf2/4/6).
  NB gsq=16 is STRONGER than any massive run (which topped at gsq=12) -> these may be marginal; WATCH
  acceptance on first jobs and bump if it drops (do not preempt).
- Packing: **2 MPS clients/GPU** (NPARALLEL_DUPDATE=1 already set for MPS).
- Accounts: **split qed3 + affine** (propose gsq=12 -> qed3, gsq=16 -> affine; each account then runs
  L2+L4 across Nf2/4/6 = balanced).
- Drivers: the **_fermilab** copies (hardcoded geometry) with **m=0.0** and **gsq as a CLI arg** (one
  L2 driver serves gsq=12 AND 16; likewise L4).

## Prep: heavy L4 cap 320 -> 120 (separate NM ask)

`hmc_fermilab_wmass_L4_heavy_claude.cu:239` `const int kmax=320;` -> `120` (comment 320; a 120 line
already exists in history at :235). Rebuild that binary. Heavy L4 currently @93-109, so 120 is a near
finish target; it is being killed now (below) and can resume to 120 later.

## Files (all new; `_claude`)

### Drivers (2 new; gsq is a CLI arg so each serves both couplings)
- `hmc_fermilab_wmass_L2_massless_claude.cu` = copy of `hmc_fermilab_wmass_L2_gsq12_claude.cu`
  (already n=21/k=0.001, nsteps=8, fermilab geometry) with **kmax 320 -> 620**. m via CLI (pass 0.0).
- `hmc_fermilab_wmass_L4_massless_claude.cu` = copy of `hmc_fermilab_wmass_L4_claude.cu`
  (n=21/k=0.001, nsteps 8/8/10) with **kmax 120 -> 320**. m via CLI (pass 0.0).
  (Both already build the diagonal mass from the CLI m; m=0 -> massless overlap.)

### Run dirs (checkpoints + logs; dir3 encodes gsq/mRe/L so no collision)
- gsq=12 -> `/lustre2/qed3/massless/`     (qed3 account)
- gsq=16 -> `/lustre2/affine/massless/`   (affine account)

### Scripts (mirror the HMC MPS submission: #SBATCH + MPS daemon + 2 clients + wait + trap)
- `run_massless_mps2_claude.sh` -- 1-GPU sbatch, MPS daemon, **2 clients**; via --export: BIN(s)/L/Nf/
  gsq/OUTDIR/WALL_SEC. Each client = one ensemble (its own dir3), conn... no -- this is HMC generation,
  each client runs `hmc_..._massless.o gsq Nf nu0 0.0 0.0 MAX_SEC`. Pair (L2,L4) same Nf same gsq on
  one GPU (L2 fast + L4 slow -> balanced). Graceful wall-stop via MAX_SEC (6th arg).
- `run_wrapper_massless_claude.sh` -- per account (gsq): loop Nf{2,4,6}, submit one 2-client job per
  Nf (pairs its L2+L4), afterany-chained; jobname `hmc_ml_gsq<g>_Nf<n>`.
- `~/launch_massless_claude.sh` -- HANDOFF (dry-run default, `--smoke` 1-traj, `--apply [NCHAIN]`):
  build the 2 massless .o (+ rebuild heavy L4 @120), make run dirs, submit both accounts.

## Kill scope (before launch)

Currently 41 jobs: 36 gsq12 + 5 heavy_L4. Campaign L2/L4 already drained; no jj jobs. `scancel` all
41 (job control, not file deletion -- on-disk checkpoints kept, all resumable later). The gsq12
nsteps=8 fix + rollback persist on disk.

## Chunks (get go-ahead per chunk)

1. **Prep+kill**: edit heavy L4 cap 320->120; `scancel` the 41 other jobs. (CONFIRM before killing.)
2. **Drivers**: create the 2 massless `_fermilab` drivers (cap/nsteps/n as above); build via handoff.
3. **Scripts**: run_massless_mps2 + per-account wrappers + launch handoff.
4. **Smoke** (1 traj/ensemble) -> validate delta/admissibility/acc -> then `--apply` full.
5. Memory + (optionally) fold into check-ensembles audit.

## Open (minor, proposing)
- Account split = by gsq (12->qed3, 16->affine). Alt: by L. Using by-gsq unless told otherwise.
- Walltime per MPS job (propose 6h) + chain depth NCHAIN.
