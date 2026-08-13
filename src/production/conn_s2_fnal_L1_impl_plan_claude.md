# FNAL L=1 stride-2 connected $Y_{lm}$ pilot -- implementation plan

**Status: LAUNCHED (2026-08-10). Chunks 1-4 DONE.** Slice of
`conn_stride2_three_site_impl_plan_claude.md` (the 3-site plan). FNAL / L=1 only, on **2 qed3 GPUs**.

**RUN STATE (2026-08-10):**
- Chunk 1 (driver `_fnal` abs-geometry copy + sm_80 build) DONE -- `conn_s2/jj_conn_s2_L1_fnal_claude.o`.
- Chunk 2 (run + wrapper scripts) DONE. Two bugs found+fixed in the run script: (a) `local a=$1 b=$((a/2))`
  in one line -> `set -u` unbound-variable (split the declarations); (b) missing `source ~/env.sh`
  -> `libhdf5.so.310` not found at runtime (added, same as the HMC job).
- Chunk 3 SMOKE PASS: 16 h5 written (4 offsets x 4 cfg), conn ylm **vector+axial**+scalar l<=3 per-m,
  path `data_<esnid>/corr_ylm_conn_t00_nhits1_s1/corr.<k>.h0.h5`. **CALIBRATION: A100 L1 = ~105 s/cfg
  = SAME as TITAN V** (small L1 is latency/overhead-bound, no A100 speedup) -> L1 wall ~= **170 worker-h
  /offset ~ 7 days** on 2 A100 (4 offset-workers), NOT the 3-4 d that a 2x A100 assumption gave.
- Chunk 4 LAUNCHED: full L1 chain-8 on qed3 (`bash run_wrapper_conn_s2_fnal_claude.sh 8`). ~40 jobs
  total (~7 d) -> repeat the wrapper to top up (resumable via complete-gate); QOS=opp when qed3 drains.
  NEXT: pull-back to LOCAL (`data_*/corr_ylm_conn_t00_nhits1_s1/`) + completeness check; then L2/L3/L4
  (their own pilots; A100 may help more on the bigger lattices -- recalibrate).

## 1. Goal
Fill the stride-2 connected local-current $Y_{lm}$ tower for the **12 L=1 massless ensembles**, in place on
FNAL (`/lustre2/affine/redo`), by running the 4 disjoint residue-class offsets ($k\equiv 3,5,7,9 \bmod 10$)
as 4 MPS workers on 2 qed3 A100 GPUs. Conn only (axial $\ell=3$ is connected). No new physics -- pure CLI
`--kmin/--stride/--kmax` split of the existing driver.

Why L=1 first: cheapest ($\approx$682 worker-h $\times$ TITAN-V, all 12 done@4000/2000 so **frozen** -- no
re-run tails), and it validates the whole FNAL measurement path before the heavy L2/L3.

## 2. Key facts established (2026-08-10, on the FNAL cluster)
- **Geometry: FNAL has the full n1..n4 set** at `/project/qed3/qed3/geometry/data/` (incl n3). The base
  driver's relative path `../../geometry/data/` (`:98`) resolves ONLY from `src/production`. To run from a
  `/lustre2` workdir (output space) we use an **absolute-geometry `_fnal` copy** (1-line change, mirrors
  `hmc_hasenbusch_block_fermilab_claude.cu`).
- **Configs: all 12 L1 dirs are on `/lustre2/affine/redo`** -> measured IN PLACE, zero shipping.
- Driver = `jj_local_ylm_scalar_conn_stoch_claude.cu`: `-DN_REFINE_CLI=<L>` compile-time; args
  `--gsq --Nf --nu0 1.0 --ens-dir <dir>/ --kmin --stride --kmax --nhits 1 --t0 0 --spin-dilution`;
  massless => mass not passed (mass_re=0). Output `data_<esnid>/corr_ylm_conn_t00_nhits1_s1/corr.<k>.h<h>.h5`
  relative to CWD; `esnid = <ens_dir basename>_vmRe0.000000vmIm0.000000`. **Resumable**: per-config `complete`
  gate -> re-run skips done configs.
- LOCAL runner `run_conn_ext_claude.sh` is the shape to mirror (4 workers, 2/GPU MPS, per-L binary, MPS
  daemon mgmt) -- but its worker split is a config-RANGE split; ours is the 4 **offset** residue classes.

## 3. The 12 L1 units (one per ensemble; each runs all 4 offsets)
at=0.2 (KMAX 4000, first=1, cfg/offset=400): Nf2 g0.5, Nf4 g0.5, Nf6 g0.5, Nf2 g1.0, Nf4 g1.0, Nf6 g1.0,
Nf2 g1.5, Nf4 g1.5, Nf6 g1.5.
at=0.1 half-a_t (KMAX 2000, cfg/offset=200): Nf2 g1.0, Nf4 g1.0, Nf6 g1.0.
Offsets: `--kmin $((1+off)) --stride 10 --kmax <last+1>` for off in {2,4,6,8}. (L1 all done -> last fixed.)

## 4. Files
### Create
| file | purpose |
|---|---|
| `jj_local_ylm_scalar_conn_stoch_fnal_claude.cu` | base driver, ONLY change `:98` dir -> `/project/qed3/qed3/geometry/data/` (absolute), so it runs from a `/lustre2` workdir. |
| `run_conn_s2_fnal_claude.sh` | the sbatch body: 2 GPUs, MPS daemon up, **4 workers = 4 offsets**, loops the L1 ensembles, `CUDA_VISIBLE_DEVICES` maps 2 workers/GPU. Resumable, NO rm. Writes `data_*` into the workdir. `WORKDIR=/lustre2/affine/redo/conn_s2`. |
| `run_wrapper_conn_s2_fnal_claude.sh` | submits `run_conn_s2_fnal` to **qed3** (`--account=qed3.lq2_gpu`, 2 GPUs, 4h), afterany-chains N jobs (resume via complete-gate), like the HMC wrapper. Env QOS/WALL/NCHAIN/ACCT. |
| `tmp_build_conn_s2_L1_fnal_claude.sh` | build handoff: `-DN_REFINE_CLI=1 -arch=sm_80`, FNAL HDF5/GSL/Eigen/HighFive include+lib, into `$WORKDIR/jj_conn_s2_L1_fnal_claude.o`. |

### Modify: none (base driver untouched; split rides on CLI flags).
### End: append a Log line to `both_3d/jj_sync_blackboard_claude.md`; the LOCAL agent pulls h5 back.

## 5. Ordered chunks
- **Chunk 1 -- driver copy + build.** Make the `_fnal` copy (absolute geometry), write
  `tmp_build_conn_s2_L1_fnal_claude.sh`, NM runs it -> `jj_conn_s2_L1_fnal_claude.o` (sm_80). Verify exit 0.
  *Files:* `jj_local_ylm_scalar_conn_stoch_fnal_claude.cu`, `tmp_build_conn_s2_L1_fnal_claude.sh`
- **Chunk 2 -- run + wrapper scripts.** `run_conn_s2_fnal_claude.sh` (4 offsets x MPS 2/GPU) +
  `run_wrapper_conn_s2_fnal_claude.sh` (qed3 sbatch, chain).
  *Files:* those two.
- **Chunk 3 -- SMOKE (qos=test, 30 min).** One ensemble (Nf2 g0.5), few configs per offset, verify: h5
  written at `data_<esnid>/corr_ylm_conn_t00_nhits1_s1/corr.<k>.h<h>.h5`, `complete` present, no CUDA error,
  per-config wall printed -> real A100 s/cfg (calibrates section-2 estimate). NM runs; I read the log.
- **Chunk 4 -- launch full L1 on qed3 (2 GPU, chained), then pull-back handoff to LOCAL.**

## 6. Open questions for NM
1. **qed3 measurement GPUs**: 2 GPUs confirmed. Same `--account=qed3.lq2_gpu`, qos normal, `--time 4h`,
   MPS 2/GPU as the HMC? (I'll assume yes -- matches the HMC MPS jobs.)
2. **Workdir**: `WORKDIR=/lustre2/affine/redo/conn_s2` (output `data_*` + binary live there; configs read
   from `/lustre2/affine/redo/<ens>/`). OK, or prefer elsewhere on /lustre2?
3. Build only **L1** now (sm_80) -- L2/L3 binaries later per their own pilots. Agreed?
4. Chain depth: L1/offset ~85 A100-h est -> ~21 four-hour jobs/offset; run 4 offsets concurrently on the 2
   GPUs. Start with NCHAIN=8 and top up? (Same cadence as the HMC.)
