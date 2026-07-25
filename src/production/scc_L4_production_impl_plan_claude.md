# SCC L=4 production -- implementation plan (2026-07-16, NM)

Generate the **L=4** overlap-fermion QED3 ensembles on the **BU SCC** (`/projectnb/qfe/nmatsum/qed3`),
using the finalized L4 tuning (`src/tuning/`, shipped 2026-07-16) and the `src/production/` driver.
This is the SCC analogue of the existing local (`/mnt/barracuda22`) production launchers.

> ## CURRENT STATE (2026-07-16, latest -- read this first)
> - **Focus = Nf=2 only.** A first full submission (both sets, all Nf=2/4/6, both arches) was run and then
>   KILLED by NM (Nf4/6 too expensive; refocus on Nf2). Wrapper now defaults `ML_NF`/`MV_NF` = `2`.
> - **Per-set pairing.** MPS pairs are formed WITHIN a set, so a massless (k400) stream is never packed with
>   a massive (k60) stream. Default run (`WHICH=both`, both arches) = 7 Nf2 ensembles -> **4 GPU chains**
>   (2 massless g{2,4,6} k400 + 2 massive g6 m{0.1,0.5,1,1.5} k60), 2 MPS streams/GPU.
> - **`s_tot` (tau) 0.8 -> 1.0** for L4 (`includes/hasenbusch_ladder_claude.h:56`; old 0.8 commented above).
>   Milder Nf2 force -> longer trajectory decorrelates better. Steps {4,4,4} + gauge MG100 unchanged.
>   => **binaries must be REBUILT** (don't pass NOBUILD); the tau=0.8 binary is stale.
> - **Environment UNIFIED into `env.sh`** (repo root): `cuda/12.8` + `gcc/13.2.0` + exports `QED3_ROOT` and
>   `QED3_INC="-I$QED3_ROOT/qfe_mod/include"` (repo Eigen, installed by `set.sh`). HMC build needs ONLY
>   Eigen+CUDA; HDF5/HighFive/GSL dropped (measurement-only). `cuda/12.8` is PINNED because `cuda/13.2`
>   dropped `sm_70` (V100); 12.8 builds both sm_70 and sm_80. The batch script sources `env.sh` at RUNTIME
>   too (dynamic CUDA/gomp libs).
> - Scope knobs: `GSQ_LIST=2.0` for gsq=2 only; `WHICH=massless|massive`; `NF_LIST`, `H_RT_FIRST`(6h)/`H_RT`(12h),
>   `N_CHAIN`(4), `PE_OMP`(16), `SUBMIT_ARCHS`. Monitor: `bash tmp_monitor_L4_scc_claude.sh` (read-only).

## Sources / provenance (later timestamp wins, per NM)
- Driver base = `src/production/hmc_hasenbusch_block_claude.cu` (commit `954943f` "L4", 2026-07-16 18:18 --
  the latest block driver; the `src/both_3d` copy is older, 2026-07-15). Nf-packed block driver, `-DLREF`, mass-generic.
- Tuning header `src/tuning/includes/hasenbusch_ladder_claude.h` and `src/production/includes/...` are
  functionally IDENTICAL (same commit; only comments differ, tuning slightly newer). Production `includes/`
  is used for the build.

## Finalized L=4 tuning (from `includes/hasenbusch_ladder_claude.h`, `frozen_window_claude.h`)
- Hasenbusch ladder (M_mass coeffs): `{0.0, 0.4, 1.0}` (frame0 = physical mass; 3-frame).
- MD steps `{4, 4, 4}`, tau (`s_tot`) = **0.8**, gauge substeps **MG = 100** (gauge force dominates -> fine
  substep + short tau give the margin, incl. Nf=6; gauge substeps carry no CG so they are cheap).
- Zolotarev: frozen window (lambda_min, lambda_max) = **(0.008, 5.0)**, `n_action = 31`, `n_force = 11`
  (two-operator split-pole force; exact by Metropolis).
- L4 built with **`-DMIXED_FORCE`** (force-only mixed precision, ~1.35x; exact by Metropolis).

## SCC environment -- UNIFIED into `env.sh` (single source; source it, nothing else)
- `env.sh` (repo root) loads `cuda/12.8`, `gcc/13.2.0` and exports `QED3_ROOT` + `QED3_INC="-I$QED3_ROOT/qfe_mod/include"`.
- **The HMC/production build needs ONLY Eigen + CUDA.** Eigen is installed IN THE REPO by `set.sh`
  (`qfe_mod/include/Eigen` -> `eigen_src/Eigen`); `geodesic.h` pulls only Eigen; checkpoints are plain
  `std::ofstream`. NO HDF5 / HighFive / GSL anywhere in the driver or `includes/` (verified by grep) --
  those are MEASUREMENT (glue/jj) deps and are deliberately OUT OF SCOPE here (per NM 2026-07-16).
  (Aside: HighFive at `$HOME/opt/highfive` was built against `/share/pkg.8/hdf5/1.10.10`, so it IS consistent
  with the `hdf5/1.10.10` module -- but irrelevant to this build.)
- The FIRST build attempt failed with `Eigen/Dense: No such file` because the wrapper used the SCC module
  var `$SCC_EIGEN_INCLUDE` (empty -- the `eigen/gsl/hdf5` modules were not loading). FIX: use the repo Eigen
  via `$QED3_INC` and drop the module loads + hdf5/gsl link entirely.
- GPU arch: SCC V100 = **sm_70** (`gpu_c 7.0`), A100 = **sm_80** (`gpu_c 8.0`); build both.
- Batch: SCC SGE (`qsub`), template `src/both_3d/run_nf.sh`.

## Machine-specific changes needed (why an `_scc` driver copy)
The only per-machine edits are the GEOMETRY PATHS, which must be **ABSOLUTE** (per NM; the `.dat` files do
not exist yet -- NM will generate them once the setup is settled). Hence a driver copy, `_scc` before `_claude`
(mirrors `hmc_fermilab_claude.cu`):
- `hmc_hasenbusch_block_scc_claude.cu`:
  - line 73: `const std::string dir = "/projectnb/qfe/nmatsum/qed3/geometry/data/";` (was `"../../geometry/data/"`)
  - line 74: `#include "/projectnb/qfe/nmatsum/qed3/geometry/geodesic.h"` (was `"../../geometry/geodesic.h"`)
  - everything else identical (Comp params, `NPARALLEL_GAUGE/SORT=16`, ladder/window via headers).

## MPS: two streams per GPU (per NM)
Pack **2 processes per GPU** via the CUDA MPS daemon (`nvidia-cuda-mps-control -d`), each a separate ensemble
process sharing one GPU. Layout depends on GPU count + how many ensembles (see OPEN QUESTIONS).

## Deliverables (all in `src/production/`) -- STATUS
1. [DONE] `hmc_hasenbusch_block_scc_claude.cu` -- SCC driver copy; ONLY change vs base = absolute geometry
   paths (`const std::string dir = "/projectnb/.../geometry/data/"`, `#include "/projectnb/.../geodesic.h"`).
2. [DONE] `run_L4_scc_claude.sh` -- SGE BATCH script: one job = one GPU = two MPS-packed streams. Starts a
   private MPS daemon (per-job pipe/log dirs under `$TMPDIR`), splits `$NSLOTS` OMP threads across the 2
   streams, runs slot A (+optional slot B) via `qsub -v` params, `tee`s per-ensemble logs. Resource requests
   (`-l gpus/gpu_c/h_rt -pe omp`) come from the wrapper's qsub command line.
3. [DONE] `run_wrapper_L4_scc_claude.sh` -- LOGIN-node wrapper: module-loads toolchain, builds one binary
   per Nf per arch (`hmc_L4_scc_Nf<Nf>_<arch>.out`, `-DLREF=4 -DNF -DKMAX -DKRNG=1 -DMIXED_FORCE`), enumerates
   Nf x gsq x mass, pairs ensembles 2/GPU, qsubs one batch job per pair. `DRYRUN=1` prints qsub cmds; `NOBUILD=1`
   skips the build. Editable lists: `NF_LIST`, `GSQ_LIST`, `MASS_LIST`, `KMAX`, `ARCH`/`GPU_C`, `H_RT`, `PE_OMP`.

Build line (resolved, SCC module env vars): `-arch=$ARCH -std=c++20 -lcublas -lcusolver -lcusparse -lgomp
-Xcompiler -fopenmp -I./includes/ -I$HOME/opt/highfive/include -I$SCC_EIGEN_INCLUDE -I$SCC_HDF5_INCLUDE
-I$SCC_GSL_INCLUDE ... -L$SCC_HDF5_LIB -L$SCC_GSL_LIB -lhdf5 -lgsl -lgslcblas -lm`.

## RESOLVED (all decisions in)
- Launcher style = SGE `qsub` (batch + wrapper), NOT interactive nohup.
- L4 params confirmed by NM: 3-stage {4,4,4}, tau(s_tot)=0.8 (matches the shipped tuning header). `-DMIXED_FORCE` on.
- Provenance: later timestamp wins -> `src/production` block driver (commit 954943f) is the base.
- **GPU arch = BOTH** `sm_70` (V100, `gpu_c 7.0`) and `sm_80` (A100, `gpu_c 8.0`); separate executables per arch.
- **Ensembles = BOTH sets**: MASSLESS `Nf{2,4,6} x gsq{2,4,6} m0` KMAX=400 (9) + MASSIVE `Nf2 gsq6 m{0.1,0.5,1,1.5}`
  KMAX=60 (4) = 13 ensembles. KMAX is baked into the binary name (`_k400_`/`_k60_`) so a shared Nf gets both.
- **13 ensembles round-robined across the two arch pools** (7 -> sm_70, 6 -> sm_80), PAIRED 2/GPU via MPS.
- **Job chaining**: each pair = a DEPENDENT CHAIN of `N_CHAIN=4` links; link0 h_rt=**6h** (per NM), links1.. h_rt=12h,
  each `-hold_jid` the previous, resuming the same ensembles' checkpoints toward KMAX. `qsub -terse` captures ids.
- **Graceful time limiter ALREADY in the driver** (`max_sec` argv[6] + `wall_timer`; breaks before a trajectory
  that would exceed budget; `k_ckpoint=1` => every-trajectory checkpoint => lossless resume). The `_scc` copy
  inherits it unchanged; the batch script passes `max_sec = h_rt - BUFFER_SEC` (900s). NOTHING added to the .cu.
- `gpu_c` = 7.0 / 8.0 (BU SCC min-compute-capability syntax); NM's old `run_nf.sh` wrote `70` -- confirm if it matters.

## TO RUN (SCC login node)
```
cd /projectnb/qfe/nmatsum/qed3/src/production
DRYRUN=1 bash run_wrapper_L4_scc_claude.sh      # inspect the qsub chain (no build/submit)
bash run_wrapper_L4_scc_claude.sh               # build both-arch binaries + submit all chains
```
Prereq: `geometry/data/*.dat` must exist at RUN time (NM generates). Build does NOT need them.
Knobs: `WHICH` (massless|massive|both), `N_CHAIN`, `H_RT_FIRST`/`H_RT`, `SUBMIT_ARCHS`, `PE_OMP`, per-set NF/GSQ/MASS/KMAX.

## Build command (SCC)
```
module load cuda/12.8 gcc/13.2.0 hdf5/1.10.10 eigen/3.4.0 gsl/2.5
nvcc -arch=sm_70 -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp \
  -I./includes/ -I$HOME/opt/highfive/include \
  -I$SCC_EIGEN_INCLUDE (module) -I$SCC_HDF5_INCLUDE (module) -I$SCC_GSL_INCLUDE (module) \
  -DLREF=4 -DNF=<Nf> -DKMAX=<target> -DKRNG=1 -DMIXED_FORCE \
  hmc_hasenbusch_block_scc_claude.cu \
  -L$SCC_HDF5_LIB -L$SCC_GSL_LIB -lhdf5 -lgsl -lgslcblas -lm -o hmc_L4_scc.out
```
(Exact module include/lib env-var names to be confirmed against `module show` on SCC.)

## OPEN QUESTIONS -- ALL RESOLVED (see CURRENT STATE at top). Kept for history:
1. Ensembles: BOTH sets, Nf=2 only (Nf4/6 dropped). Targets massless 400 / massive 60.
2. GPUs: no fixed count -- one GPU per pair-chain (7 ensembles -> 4 chains), 2 MPS streams/GPU.
3. Launcher: SGE qsub (batch `run_L4_scc_claude.sh` + wrapper `run_wrapper_L4_scc_claude.sh`).
4. `-DMIXED_FORCE`: yes, stays on.

## UPDATE 2026-07-17 -- RUN IS LIVE + healthy (chronological log of what happened)
- **Geometry blocker (fixed).** First launches crashed at startup with `std::out_of_range: map::at` inside
  `D.update`. Root cause: the spin-connection files `omega_n4.dat` + `alpha_n4.dat` (read by
  `dirac_simp.h::SpinStructureSimp`) were MISSING at run time -- `pts/links/nns/faces/duals` alone are NOT
  enough. NM generated them; verified complete (omega 480 lines, alpha 960 = both link orientations, 0
  missing). NOT a compile-flag issue (SCC flags == local: `-g -O3 -std=c++20 -lcublas -lcusolver -lcusparse
  -lgomp -Xcompiler -fopenmp`). => **geometry gen for L must emit `omega_n<R>.dat` + `alpha_n<R>.dat`.**
- **Compile flags:** stay plain (NM decided): NO `-Ofast`, NO fast-math (host or device), NO `-march=native`
  (build on login, run on compute -> SIGILL risk; and the run is GPU-bound so the gain is nil). Only genuinely
  useful safe lever would be adding native `-arch` per exact GPU (sm_89 L40S, sm_86 A40) -- not done.
- **Throughput:** ~1800 s/traj (~30 min). MPS 2-stream costs ~6%/stream -> ~2x aggregate (worth it).
  Acceptance ~100% (one 0.77), |dH|<=0.6 -- clean. `N_CHAIN=4` (~80 traj) is FINE for massive (60) but far
  short of massless (400 needs ~18 links); massless needs re-runs to extend.
- **Concurrent-writer incident (fixed).** Re-running the wrapper spawned a 2nd parallel chain for 2 ensembles
  (old+new both writing one ckpt dir). Cost = wasted allocation only; NO data corruption (all `ckpoint_lat`
  intact, uniform 657408 B; contiguous k). NM `qdel`'d the old remnants. FIX: wrapper now snapshots existing
  jobs (`qstat -r`) and `existing_tail()` `-hold_jid`s a new chain onto any live chain covering either
  ensemble (match by token `Nf<nf>g<gsq>m<mass>`) -- re-runs EXTEND, never collide. Verified by dry-run.
- **rng disk:** each rng ckpt ~0.5 GB; `KRNG=1` (running binary) keeps ALL -> currently 136 ckpts / 67.5 GB
  (disk 4.6 T free of 11 T, fine). Wrapper default now `KRNG=4` (future rebuilds thin). PENDING (NM: "work on
  rng shortly"): delete old rng's keeping latest ~3/ensemble (NM runs the `rm`; Claude never rm's).
- **Intermittent `CUDA error`** at first `cudaMalloc` on some shared GPU nodes (MPS/GPU-init, node-dependent);
  the `-hold_jid` chain self-recovers on the next link/node. Watch; if frequent, revisit MPS or 1 stream/GPU.
- **Current live jobs:** 4 serial chains (verified single-writer): `L470_Nf2g2`(g2+g6), `L480_Nf2g4`(g4),
  `L470_Nf2g6`(m0.1+m1.0), `L480_Nf2g6`(m0.5+m1.5). Tails to anchor future extends: 6738719/23/27/31.
- Monitor: `bash tmp_monitor_L4_scc_claude.sh`.
- **Acceptance (measured 2026-07-17):** 87-94% per ensemble (gsq2 87, gsq4 92, gsq6 87, m0.1 90, m0.5 91,
  m1.0 94, m1.5 91). maxdH 3-8 = early-thermalization spikes, not a problem.
- **rng cleanup done:** thinned to keep-every-4th + latest-per-ensemble via `rng_delete_list_claude.txt`
  (Claude lists, NM `xargs rm`); 136->37 ckpts, 67.5->18.3 GB, resume-safe. KRNG=4 thins future links.
- **Blackboard `redo_ensembles_claude.txt`:** L4 (SCC) rows now have exact massive hb tags (resc=0.259021*mRe)
  + refreshed snapshot + an RSYNC block (host scc1.bu.edu, base `.../src/production/`, pull only
  `ckpoint_lat.*`) so a local agent can pull configs. (Shared file w/ the separate AFFINE L1/L2 campaign.)
- **SCOPE clarified (NM):** L4 Nf=2 = SCC (this plan). L4 Nf=4,6 = handled at FNAL, NOT here.
- **Chains EXTENDED (NOBUILD=1):** re-ran the wrapper; anchoring `-hold_jid`'d new links onto live tails
  (6738719/23/27/31) -> ~6-7 links/chain, 26 jobs, verified single-writer (no collision). `H_RT_FIRST` bumped
  6h->12h (all links 12h now). Committed `git fa0eea7 "from SCC"` (run_wrapper + this plan).
- **FNAL sync (bidirectional git, branch develop):** FNAL added `-DL4_MDSTEP` to `includes/hasenbusch_ladder_claude.h`
  (default 4 -> Nf2 `{4,4,4}` UNAFFECTED; Nf4/6 uses 5), L4 Nf4/6 blackboard rows (p14/p15, "not-launched"),
  and an "RSYNC (FNAL->local)" block (lq.fnal.gov:/lustre2/affine/redo/). My `tau=1.0` survived the merge.
  `tuning/includes` copy DIVERGED (still tau=0.8, no L4_MDSTEP) -- harmless for SCC (builds against production/includes).
- **NM now moving to LOCAL + observation (measurement).** Configs pulled via the blackboard's two rsync blocks
  (SCC: `scc1.bu.edu:.../src/production/Nf2_*L4_hb*/ckpoint_lat.*`; FNAL: `lq.fnal.gov:/lustre2/affine/redo/`).
  SCC jobs keep filling autonomously via the `-hold_jid` chains.

## UPDATE 2026-07-18 -- massive done; switched to 1 stream/GPU (MPS OFF)
- **Massive (4 ens) COMPLETE** at k=59 (target 60; driver `for k<kmax` gives 59). Chains drained normally.
  Wrapper `WHICH` default -> `massless` (massive lines effectively commented out; `WHICH=both` re-includes).
- **`CUDA error` root cause = MPS on shared nodes.** Frequent `what(): CUDA error` at the first cudaMalloc on
  `ece`/`l40s` nodes: `-l gpus=1` -> GPU in EXCLUSIVE_PROCESS mode; the per-job `nvidia-cuda-mps-control -d` +
  client attach is unreliable there, so a stream's first context/malloc fails. Chains self-recovered via
  `-hold_jid` but burned many links. Exact code unrecoverable (CUDA_CHECK printf buffered, lost on `terminate()`).
- **FIX: 1 stream/GPU, no MPS.** Wrapper `NPACK=1` (default now) -> one ensemble per GPU, no pairing; batch
  `run_L4_scc_claude.sh` starts MPS only when `NPROC>=2` (solo runs the binary directly, prints "single stream
  (no MPS)"). `NPACK=2` restores MPS pairing. Massless resubmitted as 3 solo chains, `-hold_jid`-anchored onto
  the still-running old 2-stream links (6738719/6738722) -> clean cutover, verified single-writer.
- **Massless progress:** gsq2 k=64 / gsq4 k=65 / gsq6 k=62 of 400 (~16%); acc 83-90%. Needs repeated re-batches
  (`NOBUILD=1 bash run_wrapper_L4_scc_claude.sh`, anchors) to reach 400.
- **1-stream/no-MPS VERIFIED (2026-07-18):** the 3 solo massless chains run clean -- logs show "single stream
  (no MPS)", ZERO CUDA errors post-switch, ~1630-1810 s/traj (no throughput loss), progress k~69-72/400.
- **MASS ERROR fixed (2026-07-18):** L4 massive was at `{0.1,0.5,1.0,1.5}` (from the STALE params_massive doc)
  != L1/L2's `{0.1,0.2,0.3,0.4}` -> WRONG physical masses. Cause: followed the doc, didn't cross-check the live
  L1/L2 blackboard rows. Fixed in wrapper (`MV_MASS="0.2 0.3 0.4"` -- only the new ones, m=0.1 already correct/done),
  `params_massive_claude.md`, and the blackboard. Old `{0.5,1.0,1.5}` dirs = wrong-mass leftovers (NM keep/discard TBD).
  NEXT: `WHICH=massive NOBUILD=1 bash run_wrapper_L4_scc_claude.sh` generates 0.2/0.3/0.4.
  LESSON: verify physics params (masses, couplings) against the ACTUAL runs, not a params doc that may be stale.

## UPDATE 2026-07-18 (late) -- massive complete, massless ~42%
- **Massive corrected set {0.1,0.2,0.3,0.4} ALL DONE** at k=59 (acc 85-92%). hb tags (resc=0.259021*mRe):
  0.1->hb0.425902-1.025902, 0.2->hb0.451804-1.051804, 0.3->hb0.477706-1.077706, 0.4->hb0.503608-1.103608.
  The wrong-mass {0.5,1.0,1.5} dirs (k=59) remain on disk as leftovers -- NM keep/discard TBD (Claude won't rm).
- **Massless ~42% of 400:** gsq2 k=167, gsq4 k=166, gsq6 k=174 (acc 87-89%), 3 solo chains, 1 stream/GPU,
  no MPS, ZERO CUDA errors. Extend to 400 via repeated `NOBUILD=1 bash run_wrapper_L4_scc_claude.sh` (anchors).
- Blackboard `redo_ensembles_claude.txt` L4 (SCC) rows updated to match (massive DONE + exact hb tags; massless ncfg).
