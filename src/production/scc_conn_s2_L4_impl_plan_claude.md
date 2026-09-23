# SCC L4 conn stride-2 run -- impl plan (_claude, 2026-08-10)

SCC's slice of the three-site stride-2 connected $Y_{lm}$ completion
(parent: `conn_stride2_three_site_impl_plan_claude.md`). **SCC does the L4 ensembles ONLY** (per NM
2026-08-10) -- SCC has complete `n4` geometry; it LACKS `omega/alpha` for n1/n2 and all of n3, so L1/L2/L3
were reassigned to LOCAL/FNAL.

## Scope: 3 ensembles x 4 offsets = 12 units (from `conn_stride2_assign_claude.txt`, site=SCC & L4)
| ensemble dir | Nf | gsq | offsets (kmin) | live last |
|---|---|---|---|---|
| Nf2_gsq2.000000...nt128L4_hb0.400000-1.000000 | 2 | 2.0 | +2/+4/+6/+8 (kmin 3/5/7/9) | ~561 |
| Nf6_gsq4.000000...nt128L4_hb0.400000-1.000000 | 6 | 4.0 | 3/5/7/9 | ~377 |
| Nf6_gsq6.000000...nt128L4_hb0.400000-1.000000 | 6 | 6.0 | 3/5/7/9 | ~276 |
Each unit = one residue class $k \equiv \text{off}{+}1 \pmod{10}$: `--kmin <off+1> --stride 10 --kmax <last+1>`.
The stride-2 target grid is the union of the existing stride-10 (k=1,11,..) + these 4 disjoint classes.

## Driver + build (NO code change; NO _scc copy)
- Driver `jj_local_ylm_scalar_conn_stoch_claude.cu` (already in `src/production/`). `N_REFINE` is COMPILE-TIME
  `-DN_REFINE_CLI=4`. Geometry dir `../../geometry/data/` (line 98) resolves from `src/production` -> no path
  edit / no `_scc` copy (unlike the HMC driver). Builds with `-I./includes/` (same as LOCAL run_conn_ext).
- **SCC build deps:** `source env.sh` (cuda/12.8, gcc, `$QED3_INC` = repo Eigen) + `module load hdf5/1.10.10 gsl`;
  HighFive `-I/projectnb/qfe/nmatsum/opt/highfive/include`. Two arches (V100=sm_70, A100=sm_80):
  ```
  nvcc -arch=sm_{70,80} -g -O3 -std=c++20 -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp \
    -DN_REFINE_CLI=4 -I./includes/ $QED3_INC -I/projectnb/qfe/nmatsum/opt/highfive/include \
    -I$SCC_HDF5_INCLUDE -I$SCC_GSL_INCLUDE \
    -L$SCC_HDF5_LIB -L$SCC_GSL_LIB -lhdf5 -lgsl -lgslcblas -lm \
    jj_local_ylm_scalar_conn_stoch_claude.cu -o jj_conn_s2_L4_<arch>.o
  ```

## Invocation (mirrors LOCAL run_conn_ext EXACTLY -- required for cross-site merge)
```
./jj_conn_s2_L4_<arch>.o --gsq <g> --Nf <nf> --nu0 1.0 --ens-dir <ensdir>/ \
  --kmin <off+1> --stride 10 --kmax <last+1> --nhits 1 --t0 0 --spin-dilution
```
Output: `data_<ESNID>/corr_ylm_conn_t00_nhits1_s1/corr.<k>.h0.h5` (machine-independent path -> plain rsync back).
Seed is deterministic (`ESNID_k<k>_h<h>`), so SCC results are poolable with the other sites (agreement at
CG TOL 1e-8, NOT bit-identical -- validate numerically, Chunk 5 of the parent).

## Scripts (this task)
- `run_wrapper_conn_s2_scc_claude.sh` -- login-node wrapper: builds the sm_70+sm_80 L4 conn binaries, reads the
  SCC&L4 units from the assignment table, computes live `kmax=last+1` per ensemble, round-robins units across
  arches with the **FP64 `gpu_type` pin (sm_70->V100, sm_80->A100)**, and submits a dependent CHAIN per unit
  (complete-gating makes each 12h link resume where the last left off). Mirrors `run_wrapper_L4_scc_claude.sh`.
- `run_conn_s2_scc_claude.sh` -- SGE batch: one job = one GPU = one (ensemble, offset) unit. Sources env +
  `module load hdf5/1.10.10 gsl` (runtime libs), runs the driver, tees a log. No submit/build/rm.
- `tmp_build_conn_s2_scc_claude.sh` -- build smoke test (compile both arches only; NM runs it, reads the log).

## Why this is safe (measurement, not generation)
- Reads configs only -> NO fork hazard with the LIVE L4 HMC on the same 3 ensembles (unlike the 2026-07-22
  HMC handoff). The conn jobs SHARE the ~6 qfe V100 slots with the running HMC -> they queue/time-share.
- Resumable + idempotent: `complete`-gated per config (atomic `.tmp`+rename); a killed/overlapping job costs a
  stat, not a solve. Re-run the wrapper after each L4 HMC top-up to sweep the new tail (`--kmax` grows).
- FP64: same weak-L40S trap as the HMC -> the wrapper pins `gpu_type=V100/A100` (`gput_of`). No L40S.

## Decisions (NM 2026-08-10)
- **L4 only on SCC** (geometry gap on n1/n2/n3). L1/L2/L3 -> LOCAL/FNAL (reweight the assignment generator).
- nhits 1, t0 0, spin-dilution ON, stride 10 per offset -- identical to the LOCAL stride-10 run.
- PE_OMP small (4) + no MPS (1 job/GPU), same as the HMC.

## Refs
Estimator: `conserved_current_correlators_impl_plan_v3` + driver header. A. D. Kennedy hep-lat/0402038 (overlap).

## Second stochastic hit (h1) -- added 2026-08-29 (NM)
Motivation: the conn measurement is far cheaper than the HMC (bandwidth-bound on idle L40S), so once an ensemble
is at fixed statistics (N_cfg frozen), the only lever left to shrink the error is the STOCHASTIC noise -> add a
2nd hit. Whether it actually helps depends on the stochastic-vs-gauge noise ratio (spin dilution already
suppresses stochastic noise, so the l=3 axial conn may already be gauge-dominated) -- but the cost is low and
N_cfg is fixed, so NM opted to just add h1 across all 3 SCC ensembles.

Error model: sigma^2 ~ (sigma_gauge^2 + sigma_stoch^2 / N_hit) / N_cfg. A 2nd hit only halves the stochastic
piece; gauge piece is untouched. Gate the go/no-go on the l=3 observable's single-hit vs 2-hit error.

Design (separate scripts + one small, default-preserving binary flag). SAME-DIR output, chosen after LOCAL
(barracuda22, owns the SCC->local rsync) flagged that its pull glob is hardcoded to `nhits1_s1/corr.*.h5` and
does NOT preserve links / pull a hit-tagged dir. So h1 goes in the SAME nhits1 dir as h0 (the disc convention),
no separate dir, no hardlinks, no rsync change (NM 2026-08-29):
- Binary flag `--outdir-nhits <N>` added to `jj_local_ylm_scalar_conn_stoch_claude.cu` (default -1 => tag the
  output dir by the run's nhits, i.e. UNCHANGED default behavior). `dir_nhits = outdir_nhits>=0 ? outdir_nhits :
  nhits` feeds the `_nhits<..>_` path token. REQUIRES REBUILD of `jj_conn_s2_L4_sm_{70,80}.out`.
- `run_conn_s2_hit2_scc_claude.sh` (batch) runs `--nhits 2 --outdir-nhits 1` -> dir_out = existing
  `corr_ylm_conn_t00_nhits1_s1/`; the per-hit "complete"-gate finds h0 there (verified: h0 files carry a
  top-level `complete` dataset) and SKIPS it -> only h1 (`corr.<k>.h1.h5`) is solved, written alongside h0.
- `run_wrapper_conn_s2_hit2_scc_claude.sh` (wrapper): all 3 SCC ensembles x 4 offsets = 12 units, job namespace
  `c2h<arch>`, default GPUT_SM80=L40S. NO pre-linking step (h0 already in the target dir).
- h1 RNG = the EXISTING per-hit deterministic seed `esnid_k<k>_h1` (same scheme as h0's `..._h0`, hit index 1 ->
  INDEPENDENT, reproducible, poolable with FNAL/LOCAL). No RNG code change.
- LOCAL pull is UNCHANGED: its `nhits1_s1/corr.*.h5` glob now also matches `corr.<k>.h1.h5`. FNAL matches by
  passing the same `--outdir-nhits 1` when it adds its h1.
- Idempotent/resumable: re-run to extend (complete-gate skips finished (k,hit) pairs; c2h anchoring extends).

Launch (NM runs):  # first REBUILD (adds the --outdir-nhits flag):
                   #   run_wrapper_conn_s2_scc_claude.sh with FORCE_BUILD=1, or the tmp build script
                   DRYRUN=1 bash run_wrapper_conn_s2_hit2_scc_claude.sh   # preview qsub lines
                   bash run_wrapper_conn_s2_hit2_scc_claude.sh            # submit 12 units (L40S), h1 -> nhits1 dir
