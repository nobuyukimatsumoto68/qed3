# Handoff: packing multiple `jj_corr_dilute` jobs onto one A100 80G (SLURM)

## Goal

A single `jj_corr_dilute_claude.cu` process leaves the GPU only ~50-60% utilized
(per `nvidia-smi`). On an **A100 80G** there is plenty of memory headroom, so we
want to **run several processes concurrently on the one allocated GPU** to raise
aggregate throughput across a config sweep. This is process-level parallelism
(orthogonal to the in-process stream concurrency), and it is the cheap, no-code
win for this measurement.

This is a throughput optimization for a many-config sweep. It does **not** speed
up any single config; wall-clock per trajectory is unchanged.

## Why a single process underfills the GPU

The problem is small and the pipeline has serial gaps that no single-process knob
removes:

1. **Tiny kernels.** At `N_REFINE=1` the full spinor length is
   `N = NS*N_SITES*Nt = 2*12*128 = 3072` (`jj_corr_dilute_claude.cu:79-88`).
   Every kernel underfills the A100 SMs; throughput is launch-latency / occupancy
   bound, not compute bound.
2. **Per-apply join barrier.** The multishift overlap apply spreads the ~21
   Zolotarev poles across `NSTREAMS` CUDA streams
   (`includes/overlap_wmass_claude.h:352, 528, 591`) but **joins on a full
   `cudaStreamSynchronize` after every apply** (`:361, :398, :538`). The device
   drains to empty between applies.
3. **Sequential CG recurrence.** Outer CG iteration `k+1` depends on `k`, so
   successive applies cannot overlap iteration-to-iteration. The barriers stack.
4. **Sink-loop CPU serialization.** The sink passes do GPU `apply_k_block_t` ->
   synchronous D2H copy -> host inner-product loop
   (`jj_corr_dilute_claude.cu:593-601`, and the parity / sp mirrors at
   `:635-639`, `:643+`). The GPU is idle during the host `dag()` accumulation.

`NSTREAMS` (currently `4`, `jj_corr_dilute_claude.cu:75`) only parallelizes work
*within* one apply. It cannot fill the inter-apply barrier, the sequential CG
stalls, or the CPU-side sink work. A **second/third/fourth process** has work
ready precisely during those gaps -- each process has its own barriers and its
own CG recurrence, so they interleave and fill each other's idle time.

> Note on a `NSTREAMS=1` thought experiment: lowering `NSTREAMS` removes
> in-process pole concurrency, making each process underfill the GPU even more,
> which leaves *more* headroom for additional processes. Process-level packing
> then substitutes for stream-level concurrency and often overlaps better
> (no shared per-apply barrier). Either way, the recommendation below stands.

## Why concurrent processes are safe (output isolation)

Output is **one file per (config, hit)** with atomic publish and a resume
sentinel:

- Path: `data_<ESNID>/corr_dil_nt0<N>_nhits<H>[_s<S>td<TD>]/corr.<k>.h<h>.h5`
  (`jj_corr_dilute_claude.cu:20, :477`).
- Each file is written to a `.tmp` then `rename`d atomically, with a `complete`
  sentinel written last (`:936-938`).
- On startup a completed file is skipped (`:478-480`).

So as long as the processes cover **disjoint config indices `k`**, they never
touch the same file and never race. Do **not** run two processes over the same
`k` range relying on resume-skip -- they would both pick the same `k` and clash
on the `.tmp`.

## How to partition the sweep

To cover an effective stride-`S` set with `P` processes, give process
`p = 0..P-1`:

```
--stride (S*P)   --kmin (KMIN + p*S)   --kmax KMAX
```

Argument parsing supports `--kmin / --kmax (exclusive) / --stride`
(`jj_corr_dilute_claude.cu:178-198`).

Worked example for the current production sweep
(`run_jj_dilute_prod_claude.sh`: `KMIN=40`, `STRIDE=8`, 140 configs,
`KMAX=1160`), split into `P=4` processes while keeping the same stride-8 set:

| process `p` | `--stride` | `--kmin` | covers k = ...            |
|-------------|-----------:|---------:|---------------------------|
| 0           | 32         | 40       | 40, 72, 104, ...          |
| 1           | 32         | 48       | 48, 80, 112, ...          |
| 2           | 32         | 56       | 56, 88, 120, ...          |
| 3           | 32         | 64       | 64, 96, 128, ...          |

All four use `--kmax 1160`. Union = `k = 40, 48, 56, ..., 1152` (the original
140 configs), partitioned with zero overlap.

## SLURM mechanics (one GPU, several processes in one job)

1. **One GPU allocation, multiple background processes.** Request
   `--gres=gpu:1`. SLURM sets `CUDA_VISIBLE_DEVICES` to your one device; every
   child process inherits it and shares that GPU. Launch the `P` instances with
   `&` and `wait`.
2. **Do NOT override `CUDA_VISIBLE_DEVICES` per process.** From the job's view
   there is exactly one visible device (index `0`); the binary already calls
   `cudaSetDevice(0)` (`jj_corr_dilute_claude.cu:256`). Setting it to `1` for a
   second process would escape the allocation and fail.
3. **Size CPUs and host memory for `P` processes.** Each process drives
   `NSTREAMS` host threads in the pole launch (the `#pragma omp parallel for
   num_threads(nstreams)` regions force `NSTREAMS` threads regardless of
   `OMP_NUM_THREADS`). Request roughly `--cpus-per-task >= P*NSTREAMS` (plus a
   little) to avoid oversubscription. Host RAM must hold `P` copies of the host
   staging buffers.
4. **MPS (recommended for small kernels).** Start a per-job MPS daemon so the
   `P` contexts share the SMs *concurrently* rather than time-slicing. Keep the
   pipe/log dirs **job-local** (keyed on `$SLURM_JOB_ID`) so concurrent jobs on a
   shared node do not collide on the default `/tmp/nvidia-mps`. On nodes set to
   compute mode `EXCLUSIVE_PROCESS`, MPS is not just a perf nicety -- it is what
   makes multi-process sharing legal at all (the daemon owns the one context;
   clients funnel through it).
5. **GPU memory budget on the A100 80G.** At `L=1` a single process is small
   (order ~1-2 GB; pool ~189 MB + `KBlockScratch` ~19 MB + per-stream
   `d_MemorySets[NSTREAMS]` + `BlockedMat` block scratch of width
   `NSTACK_TP=12` / `NSTACK_SP=30`). **Measure the real footprint first** (run
   one process, read `nvidia-smi`), then set `P = floor(75 GB / per_proc_GB)`
   with margin. 80 GB leaves room for many processes; the practical ceiling is
   diminishing returns, not memory.

## How many processes?

- `P=2` is the clear, safe win (fills the barrier / CPU-sink gaps).
- More processes keep helping until the SMs saturate; with MPS on the A100 you
  can likely push to `P=4` or higher before returns flatten.
- Empirically tune: bring up `P`, watch `nvidia-smi` utilization climb toward
  ~90-100% and check that total throughput (configs/hour summed over processes)
  actually rises. Stop when utilization saturates or aggregate throughput stops
  improving. Going past that only adds contention and host-thread pressure.

## Sketch `sbatch` script

Match the build flags / includes / `ENS_DIR` from `run_jj_dilute_prod_claude.sh`.
Skeleton (fill in `P`, account, partition, time, mem to taste):

```bash
#!/usr/bin/env bash
#SBATCH --job-name=jjdil_mps
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16        # >= P*NSTREAMS ; tune to P
#SBATCH --mem=64G                 # P copies of host buffers
#SBATCH --time=48:00:00
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

# CUDA_VISIBLE_DEVICES is already set by SLURM to the one allotted GPU; do NOT touch it.

# ---- per-job MPS daemon (concurrent SM sharing; required on EXCLUSIVE_PROCESS nodes) ----
export CUDA_MPS_PIPE_DIRECTORY=/tmp/mps_${SLURM_JOB_ID}
export CUDA_MPS_LOG_DIRECTORY=/tmp/mpslog_${SLURM_JOB_ID}
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
nvidia-cuda-mps-control -d

# ---- sweep partition: P processes over the stride-8, KMIN=40 set ----
ENS_DIR="Nf2_gsq2.000000at0.200000nu01.000000nt128L1/"
KMIN=40; S=8; KMAX=1160; P=4
COMMON="--ens-dir $ENS_DIR --kmax $KMAX --nhits 1 --n-t0 2 --time-dilution 2 --spin-dilution"

for p in $(seq 0 $((P-1))); do
  kmin_p=$(( KMIN + p*S ))
  stride_p=$(( S*P ))
  ./jj_corr_dilute.o $COMMON --kmin "$kmin_p" --stride "$stride_p" \
      > jj_dilute_mps_p${p}_claude.log 2>&1 &
done
wait

echo quit | nvidia-cuda-mps-control   # tear down the daemon
```

Notes:
- Build the binary once (same `nvcc` line as `run_jj_dilute_prod_claude.sh`)
  before the parallel launch, or have the script compile first and exit on
  failure.
- One log per process (`jj_dilute_mps_p<p>_claude.log`) so progress is readable
  per stream.
- No `rm` anywhere. If a `complete`-gated output file blocks a rerun, surface it
  and let the human delete it -- do not script the delete.

## Validation checklist

1. **Correctness unchanged.** Concurrent processes write disjoint files; the
   physics per config is identical to the serial run. Spot-check one `k` against
   a serial-run `corr.<k>.h0.h5` (bitwise-equal data under `h0/`).
2. **No file collisions.** Confirm each process's `--kmin/--stride` set is
   disjoint (table above). Two processes must never share a `k`.
3. **Utilization rose.** `nvidia-smi` during the run should sit well above the
   ~50-60% single-process baseline (target ~90%+ with MPS).
4. **Throughput rose.** Sum configs-completed/hour across processes; it should
   scale up sub-linearly with `P` and flatten once the SMs saturate.
5. **Memory safe.** Peak GPU memory (all processes) stays under 80 GB with
   margin; peak host RAM under the `--mem` request.

## Key file references

- `jj_corr_dilute_claude.cu:75` -- `NSTREAMS` (intra-apply pole concurrency).
- `jj_corr_dilute_claude.cu:256` -- `cudaSetDevice(0)` (uses the one visible GPU).
- `jj_corr_dilute_claude.cu:178-198` -- `--kmin / --kmax / --stride` parsing.
- `jj_corr_dilute_claude.cu:477-480, :936-938` -- per-(k,h) atomic file + resume
  skip (makes disjoint partitions safe).
- `jj_corr_dilute_claude.cu:593-601, :635-639` -- sink-loop GPU/CPU serialization
  (the idle gaps a second process fills).
- `includes/overlap_wmass_claude.h:352, 361, 528, 538` -- pole streams + per-apply
  sync barrier.
- `run_jj_dilute_prod_claude.sh` -- canonical build flags, includes, `ENS_DIR`,
  and the production sweep this handoff parallelizes.
