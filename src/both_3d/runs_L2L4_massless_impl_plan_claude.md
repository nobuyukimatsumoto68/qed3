# L=2,4 massless HMC SEA runs -- setup plan

## Goal

Generate massless ($m=0$) overlap-fermion HMC SEA ensembles at refinement
$L = N_\text{REFINE} = 2$ and $4$, extending the existing $L=1$ production set.
Physics params match the $L=1$ un-suffixed production ensembles:
`gsq=8.0`, `Nt=128`, `at=0.2`, `nu0=1.0`, pole count $n=11$.

## Resources / packing

Local machine, 2x TITAN V (sm_70), 12 GB each. Pack two HMC procs per GPU
co-resident under CUDA MPS (SM-shared, ~1.6-1.9x aggregate; mirrors the jj/ylm
packing convention -- MPS daemon assumed up via `nvidia-cuda-mps-control -d`).
Each binary runs `kmax=4000` trajectories then exits, resume-safe (reads
existing `ckpoint_lat.k` / `ckpoint_rng.k`).

- GPU 0: `hmc_L2_claude.o` -> Nf=2 AND Nf=4 (concurrent under MPS)
- GPU 1: `hmc_L4_claude.o` -> Nf=2 AND Nf=4 (concurrent under MPS)
- Nf=6: deferred (add `6` to `NF_LIST`).

## Files

- `hmc_L2_claude.cu` -- copy of `hmc_claude.cu` with `N_REFINE=2` (line 65). Nothing else changed.
- `hmc_L4_claude.cu` -- copy of `hmc_claude.cu` with `N_REFINE=4` (line 65).
- `run_nf_L2L4_massless_claude.sh` -- builds both binaries (sources `../../env.sh`),
  then launches the two GPU streams in parallel, each running Nf=2 then Nf=4
  sequentially; tees per-job logs.

## Notes / prerequisites (verified)

- Geometry data for `n2` and `n4` already present in `../../geometry/data/`.
- Output dir naming carries `L<N_REFINE>` (`hmc_claude.cu:255`), so L2/L4 dirs
  do NOT collide with the existing `...L1` ensembles.
- Massless: `mass = (0,0)` (`hmc_claude.cu:217`); $n=11$ poles (`hmc_claude.cu:219`).
- Build: local `Makefile` auto-globs `*.cu` (excluding `_fermilab`); targets
  `hmc_L2_claude.o`, `hmc_L4_claude.o`. Flags `-arch=sm_70 -O3 -std=c++20`.

## Open questions

(none -- packing, params, pole count resolved with user 2026-06-19)
