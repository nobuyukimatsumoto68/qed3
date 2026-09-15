# conn_s2_local_impl_plan_claude.md

Local (barracuda22) stride-2 connected $Y_{lm}$ tower generation.
Created 2026-08-28. Blackboard: `redo_ensembles_claude.txt` UPDATE 2026-08-28.

## Goal / physics

Fill the **stride-2 connected $Y_{lm}$ scalar+vector current-correlator tower** on the massless
ensembles, generated **locally on barracuda22 GPU0, sequential**, in this order:

1. **L3 at=0.2** (9 ensembles, `Nf{2,4,6}` x `gsq{1.5,3.0,4.5}`, configs local @999)
2. **L1 at=0.1** (9 ensembles, `Nf{2,4,6}` x `gsq{0.5,1.0,1.5}`, configs local @1999)

"Stride-2" = the odd-$k$ grid $k \equiv 1,3,5,7,9 \pmod{10}$. The existing conn tower already holds
$k \equiv 1 \pmod{10}$ (made by `run_conn_ext_claude.sh`, stride-10 anchored at the first config).
This job adds the **four new residue classes** $k \equiv 3,5,7,9 \pmod{10}$ (driver `--kmin 3|5|7|9
--stride 10`), whose union with the existing $k\equiv1$ class is the full stride-2 (odd-$k$) grid.

## Algorithm source / provenance

- Driver: `jj_local_ylm_scalar_conn_stoch_claude.cu` (the LOCAL, relative-geometry driver; reads
  `../../geometry/data/`). NOT the `_fnal` copy, which hardcodes `/project/qed3/qed3/geometry/data/`
  and cannot run on barracuda22. Both drivers carry the 2026-08-21 valence-$a_t$ fix (`at_from_ensdir`
  + `--at` override).
- Residue-class (stride-2) scheme mirrors `run_conn_s2_fnal_claude.sh` (offsets `2 4 6 8` -> kmin
  `3 5 7 9`); here run locally, single-GPU, sequential rather than SLURM+MPS.
- Stochastic $Z_2$ sources with spin dilution, `--t0 0`, `--nhits 1` (matches the existing conn tower).
- Block-averaging convention: source-origin blocks, arXiv:0804.1501 (as elsewhere in this project).

## Current status (live survey 2026-08-28; conn = have / stride-2 target ~ncfg/2)

L3 at=0.2 (~500 target): Nf2 g1.5=203 g3.0=196 g4.5=109 ; Nf4 g1.5=33 g3.0=130 g4.5=34 ;
  Nf6 g1.5=130 g3.0=48 g4.5=0(NOT STARTED). Short ~300-500 each.
L1 at=0.1 (~1000 target): Nf2 g0.5=400 g1.0=571 g1.5=600 ; Nf4 g0.5=400 g1.0=800 g1.5=547 ;
  Nf6 g0.5=516 g1.0=400 g1.5=274. Short ~200-726 each.

Configs all LOCAL (L3 @999, L1 at=0.1 @1999). Geometry n3 (13 families) + n1 (23) present.
Build recipe verified from `run_conn_ext_claude.sh` (local nvcc + HDF5/Eigen/HighFive include paths).

## Files to create / modify

- CREATE `run_conn_s2_local_claude.sh` -- non-SLURM, GPU0, sequential stride-2 conn runner.
- REUSE binary `jj_local_ylm_scalar_conn_stoch_L{3,1}.o` (same name/recipe as `run_conn_ext`; built
  by the script's `need_build` gate if absent or stale).
- After the run: update `redo_ensembles_claude.txt` (the UPDATE 2026-08-28 table) + this plan.

## Ordered chunks

### Chunk 1 -- build the local per-L conn binaries
Files: `run_conn_s2_local_claude.sh` (build section, copied from `run_conn_ext_claude.sh:30-71`).
Build `jj_local_ylm_scalar_conn_stoch_L3.o` and `_L1.o` with `-DN_REFINE_CLI=L`, `need_build` gate.

### Chunk 2 -- write the sequential stride-2 runner
Files: `run_conn_s2_local_claude.sh`.
- `OMP_NUM_THREADS=4`, MPS 2-pack (NPACK=2) on GPU0 -> two driver processes packed under MPS (per NM).
- OFFSETS: L3 = `1 3 5 7 9` ; L1 at=0.1 = `1 3 5 7 9` (FNAL is OUT of conn -> local owns the FULL
  odd-k grid; a per-residue survey found k=1 mod10 count = 0 on ALL 9 L1 at=0.1 ensembles).
- PHASE A: L3 at=0.2 ensembles (glob `Nf*at0.200000*nt128L3_hb*`, exclude `_vmRe`), pass `--at 0.200000`.
- PHASE B: L1 at=0.1 ensembles (glob `Nf*at0.100000*nt128L1_hb*`), pass `--at 0.100000`.
- Per ensemble, per offset: `CUDA_VISIBLE_DEVICES=0 ./jj_..._L${L}.o --gsq --Nf --nu0 1.0 --at ${AT}
  --ens-dir DIR/ --kmin ${off} --stride 10 --kmax $((last+1)) --nhits 1 --t0 0 --spin-dilution`.
- Resumable: driver per-config "complete" gate skips existing h5. NO rm anywhere.
- Per-ensemble/offset logs `conn_s2L${L}_g..._Nf..._off..._claude.log`; top log via nohup tee.

### Chunk 3 -- handoff
Files: `run_conn_s2_local_claude.sh` (self-contained: builds then runs).
User runs detached:
  `nohup bash run_conn_s2_local_claude.sh > run_conn_s2_local_claude.log 2>&1 &`
(GPU MEASUREMENT job -> user launches, per the run-approval convention.)

### Chunk 3b -- glueball REMAINDER (CPU, concurrent, NO new code)
Files: (none new) re-run existing `run_glue_f2_v2_sweep_claude.sh`.
The v2 sweep loops L1-4, selects all massless ensembles at stride 1, and is complete-gated -> a re-run
picks up EXACTLY the leftover: L1 at=0.1 (Nf4/Nf6 g0.5/g1.5, un-run) + L4 partial. CPU-only
(`OMP_NUM_THREADS=1`, NWORK workers) -> runs concurrently with the GPU conn, no contention. Launch in
its own terminal:
  `nohup bash run_glue_f2_v2_sweep_claude.sh > run_glue_f2_v2_sweep_claude.log 2>&1 &`

### Chunk 4 -- reconcile
After completion: re-survey conn + glue counts, update the blackboard table + per-row `obs` tags.

## Resolved decisions (2026-08-28, NM)

- **k=1 topup**: YES for BOTH L3 and L1 at=0.1 (offsets `1 3 5 7 9`). FNAL is out of conn (allocation
  exhausted -> HMC-only now), so local owns the full odd-k grid; k=1 mod10 was found entirely absent (0)
  on all 9 L1 at=0.1 ensembles, so PHASE B must include offset 1 (revised from the initial `3 5 7 9`).
- **Packing / GPUs**: GPU0 only, MPS NPACK=2 = 2 workers (NM reverted 2026-08-28 from the brief 2-GPU
  trial). `GPU_LIST="0 1"` env override re-enables both TITAN V (4 workers) if wanted later.
- **Glueball remainder**: INCLUDED -- re-run `run_glue_f2_v2_sweep_claude.sh` (Chunk 3b), CPU-concurrent.
- **L4** (NM 2026-08-28): L4 CONN = SCC (pull back); L4 DISC + GLUE = LOCAL. The conn runner does NOT
  touch L4 (correct). L4 glue launched first via the sweep's new `L_ONLY=4` knob. L4 disc local = TODO.
