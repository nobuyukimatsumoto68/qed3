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
