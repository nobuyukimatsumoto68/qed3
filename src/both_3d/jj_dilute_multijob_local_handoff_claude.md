# Handoff (LOCAL box): packing multiple `jj_corr_dilute` jobs per GPU

Companion to `jj_dilute_multijob_handoff_claude.md` (the SLURM / A100 80G
version). Same idea -- a single `jj_corr_dilute_claude.cu` process only reaches
~50-60% GPU utilization, so pack several processes onto one GPU to raise
aggregate throughput across a config sweep. This file is the **local
workstation** recipe (no SLURM).

## Local hardware (measured 2026-06-13)

- **2x NVIDIA TITAN V**, indices `0` and `1`, **12 GB each**, arch **sm_70**
  (Volta), compute mode **Default**.
- "Default" compute mode means multiple processes can share a GPU **without**
  MPS (they time-slice). MPS is still worth enabling for the small-kernel
  concurrency win, but it is optional here (unlike EXCLUSIVE_PROCESS nodes).
- **12 GB is tight** -- much less headroom than the A100 80G. Measure the
  per-process footprint first and keep margin; you will fit fewer processes per
  GPU than on the A100.

## Why this works / why one process underfills

Identical to the A100 handoff -- see `jj_dilute_multijob_handoff_claude.md`
"Why a single process underfills the GPU". In brief: tiny kernels at
`N = 3072` (`jj_corr_dilute_claude.cu:79-88`), a per-apply
`cudaStreamSynchronize` barrier (`includes/overlap_wmass_claude.h:361, 398`), a
sequential CG recurrence, and the sink-loop GPU->D2H->host-`dag()`
serialization (`jj_corr_dilute_claude.cu:593-601, :635-639`). A second process
has work ready during exactly those idle gaps.

## Output isolation (why concurrent processes are safe)

One file per `(config, hit)`:
`data_<ESNID>/corr_dil_nt0<N>_nhits<H>[_s<S>td<TD>]/corr.<k>.h<h>.h5`, written
`.tmp` then atomically `rename`d with a `complete` sentinel last, and skipped on
restart if complete (`jj_corr_dilute_claude.cu:477-480, :936-938`). Disjoint
config sets => no shared files, no races. Never run two processes over the same
`k` range relying on resume-skip.

## Two ways to spread work on this box

1. **Across the two physical GPUs** (first, free win): one process on GPU 0, one
   on GPU 1, via `CUDA_VISIBLE_DEVICES`. This already ~2x's throughput with zero
   sharing.
2. **Multiple processes per GPU** (this handoff): stack `P` processes on each
   GPU to soak up the ~40-50% idle. Combine both: `G` GPUs x `P` procs/GPU.

## Partition scheme

For `M = G*P` total processes covering an effective stride-`S` set, give process
`j = 0..M-1`:

```
--stride (S*M)   --kmin (KMIN + j*S)   --kmax KMAX
```

Argument parsing: `--kmin / --kmax (exclusive) / --stride`
(`jj_corr_dilute_claude.cu:178-198`). Example, current production sweep
(`run_jj_dilute_prod_claude.sh`: `KMIN=40`, `STRIDE=8`, `KMAX=1160`, 140
configs), 2 GPUs x 2 procs = `M=4`:

| proc `j` | GPU | `--stride` | `--kmin` | covers k = ...   |
|---------:|----:|-----------:|---------:|------------------|
| 0        | 0   | 32         | 40       | 40, 72, 104, ... |
| 1        | 0   | 32         | 48       | 48, 80, 112, ... |
| 2        | 1   | 32         | 56       | 56, 88, 120, ... |
| 3        | 1   | 32         | 64       | 64, 96, 128, ... |

Union = the original 140 configs, partitioned with no overlap.

## Ready-to-run launcher script

Self-contained: builds once (early-exit on failure), optionally starts a
per-GPU MPS daemon, launches `G*P` disjoint partitions in the background,
`wait`s, then tears MPS down. Build flags / includes / `ENS_DIR` follow
`run_jj_dilute_prod_claude.sh`. **`-arch=sm_70`** is correct for the Titan V
(do NOT use `sm_80` here -- that is the A100 handoff).

Save as e.g. `run_jj_dilute_multi_claude.sh` and run with
`bash run_jj_dilute_multi_claude.sh`.

```bash
#!/usr/bin/env bash
# LOCAL multi-process driver for jj_corr_dilute_claude.cu on the 2x TITAN V box.
# Packs G GPUs x P processes/GPU; partitions the stride-S sweep with no file
# collisions. NO rm anywhere (delete blocking outputs by hand, see note below).
set -u
cd /mnt/barracuda22/qed3/qed3/src/both_3d || exit 1

module load cuda/12.8 2>/dev/null || true
module load gcc/13.2.0 2>/dev/null || true

# ===================== build (once, early-exit on failure) =====================
NVCC=nvcc
NVCCFLAGS="-arch=sm_70 -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
INCLUDES='-I./includes/ -I/projectnb/qfe/nmatsum/qed3/opt/eigen -I/opt/eigen-3.4.0/ -I/mnt/hdd_barracuda/opt/highfive/include/ -I/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/include/'
LDFLAGS='-L/mnt/hdd_barracuda/opt/myhdfstuff/hdf5-2.1.0/lib/ -L/usr/lib/ -L/usr/local/lib/ -lhdf5 -lgsl -lgslcblas -lm'

BIN=jj_corr_dilute.o
echo "### compile jj_corr_dilute_claude.cu  [$(date +%F_%H:%M:%S)] ###"
$NVCC $NVCCFLAGS $INCLUDES $LDFLAGS jj_corr_dilute_claude.cu -o "$BIN"
st=$?
if [ "$st" -ne 0 ]; then echo "### BUILD FAILED (status $st) ###"; exit 1; fi

# ===================== layout =====================
GPUS=(0 1)          # physical GPUs to use
P=2                 # processes PER GPU (raise only after checking 12 GB headroom!)
USE_MPS=1           # 1 = per-GPU MPS daemon (SM-concurrency for small kernels); 0 = plain time-slice

ENS_DIR="Nf2_gsq2.000000at0.200000nu01.000000nt128L1/"   # bare dir; holds ckpoint_lat.k
KMIN=40
S=8
KMAX=1160           # exclusive
COMMON="--ens-dir $ENS_DIR --kmax $KMAX --nhits 1 --n-t0 2 --time-dilution 2 --spin-dilution"

G=${#GPUS[@]}
M=$(( G*P ))        # total processes
echo "### launch G=$G GPUs x P=$P = M=$M processes  [$(date +%F_%H:%M:%S)] ###"

# OpenMP threads per process: the pole-launch regions force NSTREAMS(=4) threads;
# cap host BLAS so M processes do not oversubscribe the cores.
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

declare -a MPS_PIPES=()
j=0
for g in "${GPUS[@]}"; do
  export CUDA_VISIBLE_DEVICES="$g"        # both the MPS daemon and its clients see ONLY GPU g

  if [ "$USE_MPS" -eq 1 ]; then
    export CUDA_MPS_PIPE_DIRECTORY="/tmp/mps_gpu${g}_$$"
    export CUDA_MPS_LOG_DIRECTORY="/tmp/mpslog_gpu${g}_$$"
    mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
    nvidia-cuda-mps-control -d
    MPS_PIPES+=("$CUDA_MPS_PIPE_DIRECTORY")
  fi

  for p in $(seq 0 $((P-1))); do
    kmin_j=$(( KMIN + j*S ))
    stride_j=$(( S*M ))
    echo "  proc $j : GPU $g  --kmin $kmin_j --stride $stride_j  (k = $kmin_j, $((kmin_j+stride_j)), ...)"
    ./"$BIN" $COMMON --kmin "$kmin_j" --stride "$stride_j" \
        > "jj_dilute_multi_g${g}_p${p}_claude.log" 2>&1 &
    j=$(( j+1 ))
  done
done
wait
echo "### all $M processes finished  [$(date +%F_%H:%M:%S)] ###"

# ===================== tear down MPS (per pipe dir) =====================
for pipe in "${MPS_PIPES[@]}"; do
  echo quit | CUDA_MPS_PIPE_DIRECTORY="$pipe" nvidia-cuda-mps-control
done
```

## Local-specific cautions

1. **12 GB memory ceiling.** Run with `P=1` first, read per-process GPU memory
   from `nvidia-smi`, then choose `P` so `P * per_proc_GB` fits 12 GB with
   margin. At `L=1` a single process is small (order ~1-2 GB), so `P=2-3` per
   Titan V is plausible -- but verify, do not assume. An out-of-memory process
   dies mid-run.
2. **`CUDA_VISIBLE_DEVICES` must be set BEFORE the MPS daemon starts** so the
   daemon binds to the intended GPU; its clients inherit the same value. The
   script does this inside the per-GPU loop.
3. **MPS pipe/log dirs keyed per GPU + PID** (`/tmp/mps_gpu${g}_$$`) so the two
   GPUs' daemons (and other users' jobs) do not collide on the default
   `/tmp/nvidia-mps`. Each daemon is torn down via its own pipe dir at the end.
4. **CPU oversubscription.** Each process forces `NSTREAMS`(=4) host threads in
   the pole launch. With `M` processes that is `~4*M` busy threads; make sure the
   box has the cores, or lower `P`. (`OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`
   set to 4 here cap host BLAS per process.)
5. **No `rm` in the script.** If a `complete`-gated output blocks a rerun, the
   agent must flag it and ask the human to delete it -- never script the delete.

## Notes for agents working in THIS local environment

- **Ask before compiling or running.** Per project rules, do not run
  `make` / `nvcc` / the binary unilaterally; hand the script to the human (or
  use the `tmp_claude.sh` + `tee` to a `*_claude.log` handoff pattern) and read
  the log back.
- **Resource cap on agent-run commands.** Anything the agent runs *directly*
  (foreground/background via the shell) must cap CPU to <=4 cores
  (`OMP_NUM_THREADS=4` etc.). This cap does NOT apply to this launcher script,
  which the human runs and which is sized for the real workload.
- **Validation:** spot-check one `k` against a serial-run
  `corr.<k>.h0.h5` (bitwise-equal under `h0/`); confirm `nvidia-smi` climbs above
  the ~50-60% baseline; confirm partitions are disjoint; confirm peak memory
  stays under 12 GB per GPU.

## Key file references

- `jj_corr_dilute_claude.cu:75` -- `NSTREAMS` (intra-apply pole concurrency).
- `jj_corr_dilute_claude.cu:256` -- `cudaSetDevice(0)` (uses the one visible GPU).
- `jj_corr_dilute_claude.cu:178-198` -- `--kmin / --kmax / --stride` parsing.
- `jj_corr_dilute_claude.cu:477-480, :936-938` -- per-(k,h) atomic file + resume
  skip (makes disjoint partitions safe).
- `jj_corr_dilute_claude.cu:593-601, :635-639` -- sink-loop GPU/CPU serialization.
- `includes/overlap_wmass_claude.h:352, 361, 528, 538` -- pole streams + barrier.
- `run_jj_dilute_prod_claude.sh` -- canonical build flags / includes / `ENS_DIR`.
- `jj_dilute_multijob_handoff_claude.md` -- the A100 80G / SLURM companion.
