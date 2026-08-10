# Stride-2 connected $Y_{lm}$ tower, split across three sites -- implementation plan

**Status: PLAN ONLY (2026-08-10). Nothing built, nothing launched, nothing transferred.**

## 1. Goal

Take the connected local-current $Y_{lm}$ tower from the present **stride-10** config grid to
**stride 2** on all 42 massless redo ensembles, by splitting the missing solves across three sites:

| site | machine | notes |
|---|---|---|
| LOCAL | `barracuda22`, 2x TITAN V | `/mnt/barracuda22/qed3/qed3/src/production` |
| FNAL | `nmatsum@lq.fnal.gov` | configs under `/lustre2/qed3`, `/lustre2/affine` |
| SCC | BU SCC, `/projectnb/qfe/nmatsum/qed3` | SGE `qsub`, V100 (sm_70) + A100 (sm_80) |

**Physics motivation.** The driver for this is the **axial $\ell = 3$** channel, whose constant+exp fit
(`effmass_axial_l3_expfit_claude.py`) and the ratio $\Delta_{A,\ell=3}/\Delta_{A,\ell=1}$ are the
noisiest fermionic quantities in the campaign, and whose $(a/R)^2$ extrapolation (L1 1.39 -> L4 1.87
against the free value 2) is currently limited by the L3/L4 errors. Axial is **connected only** -- the
disconnected driver computes the vector current alone
(`jj_local_ylm_disc_stoch_claude.cu:313`, the `a=0..2` loop; no axial anywhere in its h5) -- so this is
a conn-only job. Disc is untouched by this plan.

**No new algorithm.** The estimator, the driver and the output format are unchanged; this is purely a
work-partitioning exercise. The estimator itself is documented in
`conserved_current_correlators_impl_plan_v3` and in the driver header
(`jj_local_ylm_scalar_conn_stoch_claude.cu:22-30`): $Z_2$ wall source at $t_0 = 0$, spin dilution into
two single-spin classes, $\phi = D_m^{-1}\eta$, $\psi^a_{lm} = D_m^{-\dagger}\Sigma^a_{lm}(t_0)\eta$,
$$
g^a_{lm}(t-t_0) = \psi^{a\dagger}_{lm}\left[\Sigma^a_{lm}(t)\,\phi\right].
$$

**Caveat carried forward, unresolved:** stride 2 gives $5\times$ the configs but $5\times$ the
statistics only if the axial $\ell=3$ correlator decorrelates in under ~2 trajectories. Nothing on disk
can measure that (the finest sampling anywhere is stride 10). See Open Question Q1.

## 2. Why the partition is clean

Three properties of the existing driver make a multi-site split safe with no code change to the physics:

1. **Per-config independence.** Each $k$ is a self-contained solve; there is no cross-config state.
2. **The k-grid is already a CLI knob.** `--kmin/--stride/--kmax`
   (`jj_local_ylm_scalar_conn_stoch_claude.cu:152-153`) is exactly how the four local workers already
   split the stride-10 run.
3. **`complete`-gating is free.** `:346,351` -- if all hits for a $k$ are done the driver `continue`s
   *without* `U.read` or `Dm.update`. An accidental overlap therefore costs a directory stat, not a solve.
4. **Output paths are machine-independent.**
   `data_<ESNID>/corr_ylm_conn_t00_nhits1_s1/corr.<k>.h<h>.h5` with
   `ESNID = <ens_dir basename> + "_vmRe<m>vmIm<m>"` (`:307`). The same ensemble produces the same path on
   every machine, so merging back is a plain `rsync`, no renaming, no index.

### 2.1 The four disjoint sub-streams

Every ensemble has `first = 1`, and the existing conn data sits on $k = 1, 11, 21, \dots$
(i.e. $k \equiv 1 \bmod 10$). The stride-2 target grid is $k = 1, 3, 5, \dots$. The difference is
therefore **exactly four residue classes mod 10**:

$$
k \equiv 3,\ 5,\ 7,\ 9 \pmod{10}
$$

expressed with the existing flags as `--kmin $((first+off)) --stride 10` for `off = 2,4,6,8`.

**Consequence: no already-solved data needs to be shipped anywhere.** Each site's sub-stream is disjoint
from the local stride-10 set by construction. (Shipping the finished h5 would be 6.3 GB out and would buy
nothing.)

The four classes are equal in cost to better than 0.5%:

| offset | flag | configs | worker-h (TITAN V) |
|---|---|---|---|
| +2 | `--kmin 3 --stride 10` | 9377 | 976 |
| +4 | `--kmin 5 --stride 10` | 9376 | 976 |
| +6 | `--kmin 7 --stride 10` | 9373 | 974 |
| +8 | `--kmin 9 --stride 10` | 9368 | 972 |
| **total** | | **37494** | **3897** |

Per-$L$ subtotals, per offset (x4 for all four):

| $L$ | worker-h per offset | x4 |
|---|---|---|
| 1 | 170 | 681 |
| 2 | 391 | 1564 |
| 3 | 172 | 690 |
| 4 | 242 | 969 |
| all | 976 | 3904 |

Timings are **measured** TITAN V numbers under 2-proc/GPU MPS packing, extracted from the
`[... s]` markers the driver prints per config in `conn_*_claude.log` (8992 configs). They will not
transfer to V100/A100 unchanged -- see Chunk 2.

### 2.2 Per-ensemble work unit table

The natural work unit is an **(ensemble, offset) pair** -- 42 x 4 = 168 units. Assign whole units to
sites; never split one unit across sites.

| ens | $a_t$ | first | last | cfg per offset | s/cfg | worker-h per offset | x4 offsets |
|---|---|---|---|---|---|---|---|
| L1 Nf2 g0.5 | 0.2 | 1 | 3999 | 400 | 110 | 12.2 | 49.0 |
| L1 Nf4 g0.5 | 0.2 | 1 | 3999 | 400 | 107 | 11.9 | 47.5 |
| L1 Nf6 g0.5 | 0.2 | 1 | 3999 | 400 | 106 | 11.8 | 47.2 |
| L1 Nf2 g1.0 | 0.1 | 1 | 1999 | 200 | 135 | 7.5 | 30.0 |
| L1 Nf2 g1.0 | 0.2 | 1 | 3999 | 400 | 135 | 15.0 | 59.9 |
| L1 Nf4 g1.0 | 0.1 | 1 | 1999 | 200 | 149 | 8.3 | 33.2 |
| L1 Nf4 g1.0 | 0.2 | 1 | 3999 | 400 | 149 | 16.6 | 66.3 |
| L1 Nf6 g1.0 | 0.1 | 1 | 1999 | 200 | 142 | 7.9 | 31.6 |
| L1 Nf6 g1.0 | 0.2 | 1 | 3999 | 400 | 142 | 15.8 | 63.3 |
| L1 Nf2 g1.5 | 0.2 | 1 | 3999 | 400 | 188 | 20.9 | 83.8 |
| L1 Nf4 g1.5 | 0.2 | 1 | 3999 | 400 | 182 | 20.3 | 81.1 |
| L1 Nf6 g1.5 | 0.2 | 1 | 3999 | 400 | 200 | 22.2 | 88.7 |
| L2 Nf2 g1.0 | 0.2 | 1 | 3999 | 400 | 202 | 22.5 | 89.9 |
| L2 Nf4 g1.0 | 0.2 | 1 | 3999 | 400 | 221 | 24.6 | 98.3 |
| L2 Nf6 g1.0 | 0.2 | 1 | 3999 | 400 | 215 | 23.9 | 95.7 |
| L2 Nf2 g2.0 | 0.1 | 1 | 1999 | 200 | 343 | 19.1 | 76.3 |
| L2 Nf2 g2.0 | 0.2 | 1 | 3999 | 400 | 343 | 38.1 | 152.5 |
| L2 Nf4 g2.0 | 0.1 | 1 | 1999 | 200 | 318 | 17.7 | 70.8 |
| L2 Nf4 g2.0 | 0.2 | 1 | 3999 | 400 | 318 | 35.4 | 141.6 |
| L2 Nf6 g2.0 | 0.1 | 1 | 1999 | 200 | 319 | 17.7 | 70.9 |
| L2 Nf6 g2.0 | 0.2 | 1 | 3999 | 400 | 319 | 35.5 | 141.8 |
| L2 Nf2 g3.0 | 0.2 | 1 | 3999 | 400 | 450 | 50.0 | 200.1 |
| L2 Nf4 g3.0 | 0.2 | 1 | 3999 | 400 | 472 | 52.5 | 210.0 |
| L2 Nf6 g3.0 | 0.2 | 1 | 3999 | 400 | 487 | 54.1 | 216.6 |
| L3 Nf2 g1.5 | 0.2 | 1 | 799 | 80 | 600 | 13.3 | 53.3 |
| L3 Nf4 g1.5 | 0.2 | 1 | 639 | 64 | 612 | 10.9 | 43.6 |
| L3 Nf6 g1.5 | 0.2 | 1 | 527 | 53 | 642 | 9.5 | 37.8 |
| L3 Nf2 g3.0 | 0.2 | 1 | 799 | 80 | 871 | 19.4 | 77.4 |
| L3 Nf4 g3.0 | 0.2 | 1 | 526 | 53 | 930 | 13.7 | 54.8 |
| L3 Nf6 g3.0 | 0.2 | 1 | 471 | 47 | 930 | 12.1 | 48.6 |
| L3 Nf2 g4.5 | 0.2 | 1 | 799 | 80 | 1901 | 42.2 | 169.0 |
| L3 Nf4 g4.5 | 0.2 | 1 | 479 | 48 | 1777 | 23.7 | 94.8 |
| L3 Nf6 g4.5 | 0.2 | 1 | 559 | 56 | 1777 | 27.6 | 110.6 |
| L4 Nf2 g2.0 | 0.2 | 1 | 599 | 60 | 1151 | 19.2 | 76.7 |
| L4 Nf4 g2.0 | 0.2 | 1 | 503 | 51 | 1075 | 15.2 | 60.9 |
| L4 Nf6 g2.0 | 0.2 | 1 | 367 | 37 | 1097 | 11.3 | 45.1 |
| L4 Nf2 g4.0 | 0.2 | 1 | 598 | 60 | 1857 | 31.0 | 123.8 |
| L4 Nf4 g4.0 | 0.2 | 1 | 427 | 43 | 1861 | 22.2 | 88.9 |
| L4 Nf6 g4.0 | 0.2 | 1 | 377 | 38 | 1860 | 19.6 | 78.5 |
| L4 Nf2 g6.0 | 0.2 | 1 | 570 | 57 | 3425 | 54.2 | 216.9 |
| L4 Nf4 g6.0 | 0.2 | 1 | 416 | 42 | 3661 | 42.7 | 170.8 |
| L4 Nf6 g6.0 | 0.2 | 1 | 276 | 28 | 3440 | 26.8 | 107.0 |

`last` is the current `ckpoint_lat` maximum. **It grows** for L3/L4 (still generating toward 800/600),
so each site's `--kmax` should be passed as `last+1` at launch time and the sub-stream re-run after
each config top-up; `complete`-gating makes the re-run cost only the new tail.

## 3. Transfer budget

With the Q2 hybrid (machine per ensemble) each site ships only its own ensembles' configs, so the
13.6 GB worst case splits roughly into thirds. Actual, from `conn_stride2_assign_claude.txt`:

| site | ensembles | worker-h | share | `ckpoint_lat` out | conn h5 back |
|---|---|---|---|---|---|
| LOCAL | 14 | 1302 | 33.4% | (already local) | 9.3 GB |
| FNAL | 14 | 1302 | 33.4% | 4.7 GB | 7.7 GB |
| SCC | 14 | 1300 | 33.3% | 4.2 GB | 9.2 GB |
| total | 42 | 3904 | | **8.9 GB out** | **26.2 GB back** |

Already-solved conn h5 (6.3 GB) is **not** shipped -- the sub-streams are disjoint from it.

Per-config `ckpoint_lat` sizes: L1 43008 B, L2 165888 B, L3 370688 B, L4 657408 B. The conn h5 is
700456 B at **every** $L$ -- its size is set by $N_t$ and the $\ell \le 3$ per-$m$ tower, not by the
lattice.

`ckpoint_rng.*` is **not** needed: this is measurement, not chain continuation.

## 4. Files

### To create (LOCAL, in `src/production/`)

| file | purpose |
|---|---|
| `conn_stride2_three_site_impl_plan_claude.md` | this plan |
| `conn_stride2_assign_claude.py` | generates the assignment table. `WEIGHTS` at the top sets the per-site capacity share (Q7); balanced greedy partition of whole ensembles. |
| `conn_stride2_assign_claude.txt` | the authoritative assignment table: one line per (ensemble, offset) -> site, with the exact `--kmin/--stride/--kmax`. Machine-readable; consumed by the run scripts, the rsync scripts and the merge check. |
| `conn_stride2_cmds_fnal_claude.txt` | FNAL deliverable per Q4: the literal per-unit command lines, no batch script. NM wraps them for the FNAL scheduler. |
| `run_conn_s2_claude.sh` | LOCAL runner. Same shape as `run_conn_ext_claude.sh` (4 workers, 2/GPU under MPS) but takes its (ensemble, offset) units from the assignment file, and runs only the ensembles assigned to LOCAL. |
| `run_conn_s2_scc_claude.sh` | SCC SGE batch script: one job = one GPU = one (ensemble, offset) unit, `--kmin/--stride/--kmax` from the assignment file, `h_rt`-sized k-ranges. |
| `run_wrapper_conn_s2_scc_claude.sh` | SCC login-node wrapper: builds the per-$L$ conn binaries for sm_70 and sm_80, enumerates SCC's assigned units, `qsub`s one job per unit. Mirrors `run_wrapper_L4_scc_claude.sh`. |
| `rsync_push_configs_s2_claude.sh` | LOCAL -> remote `ckpoint_lat.*` push, `--files-from` driven by the assignment file, resumable, no `--delete`, **no `rm`**. NM runs it (auth). |
| `rsync_pull_conn_s2_claude.sh` | remote -> LOCAL conn h5 pull, glob `data_*/corr_ylm_conn_t00_nhits1_s1/`, no `--delete`, **no `rm`**. NM runs it. |
| `check_conn_s2_claude.py` | merge/validation: per (ensemble, offset) count the h5 present vs expected, list gaps, verify `complete`, report duplicates. |

### To create (remote-only driver copies, if Chunk 2 confirms they are needed)

| file | purpose |
|---|---|
| `jj_local_ylm_scalar_conn_stoch_scc_claude.cu` | ONLY change vs the base: absolute geometry paths. `:98` `const std::string dir = "../../geometry/data/"` -> `"/projectnb/qfe/nmatsum/qed3/geometry/data/"`. Same rationale as `hmc_hasenbusch_block_scc_claude.cu`. |
| `jj_local_ylm_scalar_conn_stoch_fnal_claude.cu` | ditto for the FNAL path. |

### To modify

**None.** `jj_local_ylm_scalar_conn_stoch_claude.cu` is not touched; the split rides entirely on
existing CLI flags. If a driver copy turns out to be unnecessary (job `cd`s to `src/production`, so the
relative geometry path resolves), the copies are dropped too.

### To update at the end

| file | change |
|---|---|
| `jj_sync_blackboard_claude.md` (in `src/both_3d/`) | append a dated Log entry announcing the three-site conn split, the assignment, and the pull-back |
| `effmass_axial_l3_expfit_claude.py`, `effmass_conn_claude.py`, `effmass_plateau_tables_claude.py` | nothing structural -- they glob `corr.*.h0.h5`, so denser data is picked up automatically. Only `KMIN` (thermalization cut, currently 20) is worth re-checking. |
| memory `project_redo_obs_flow.md`, `project_jj_fermilab_sync.md` | record the split |

## 5. Ordered chunks

### Chunk 0 -- DROPPED (Q1: stride 2 now, no $\tau_\text{int}$ pilot)

### Chunk 1 -- assignment table [DONE 2026-08-10]
Freeze the 168 (ensemble, offset) units. Per Q2 the **machine is chosen per ensemble** and all four of
that ensemble's offsets go to the same site; the offsets are the intra-site jobs. Balanced greedy
partition of whole ensembles by worker-h, capacity weights from Q7 (default equal thirds).
Format: one line `ensemble L Nf gsq at offset kmin stride kmax site worker_h`, comment header.
*Files:* `conn_stride2_assign_claude.py`, `conn_stride2_assign_claude.txt`,
`conn_stride2_cmds_fnal_claude.txt`

Result at equal-thirds weights: 14 ensembles and 1300-1302 worker-h per site (balanced to 0.1%).
**Consequence to note for Chunk 2:** the balanced partition gives every site a mix of all four $L$, so
each site needs all four binaries *and* the full `n1`..`n4` geometry -- including the new `n3`, which is
verified only LOCALLY. If a site turns out to lack `n3`, either ship the geometry `.dat` files (small)
or re-run the generator with an $L$ constraint.

### Chunk 2 -- remote site survey (NM runs; I cannot ssh)
For each of FNAL and SCC, establish and record:
- which of the 42 massless ensemble dirs exist, and their `ckpoint_lat` ranges
- whether `geometry/data/` has the `n1`..`n4` sets (`pts/links/nns/faces/duals`, `omega_n*`, `alpha_n*`)
  -- **`n3` is new and is verified only LOCALLY**; the 2026-07-22 SCC note verified `n4` only
- HDF5 / HighFive / GSL availability and the include/lib paths
- GPU model + count, batch system, walltime cap, queue policy
- one **calibration** conn run: 3 configs on a known ensemble, to convert the TITAN V s/cfg column of
  section 2.2 into that site's s/cfg
  *Files:* `conn_s2_site_survey_claude.md` (I write it from NM's pasted output)

### Chunk 3 -- build at each site
Per-$L$ binaries `jj_local_ylm_scalar_conn_stoch_L{1,2,3,4}.o`, `-DN_REFINE_CLI=L`, per-site
include/lib line. On SCC build for **both** sm_70 and sm_80. Only the $L$ actually assigned to a site
need building there.
*Files:* `run_wrapper_conn_s2_scc_claude.sh`,
`jj_local_ylm_scalar_conn_stoch_{scc,fnal}_claude.cu` (only if Chunk 2 says the relative geometry path
does not resolve). FNAL build/submission is NM's per Q4.

### Chunk 4 -- push configs
Push only the `ckpoint_lat.*` of the ensembles each site was assigned, driven by
`conn_stride2_assign_claude.txt`. NM runs it (`kinit` for FNAL kerberos first; the script only warns via
`klist -s`, it never runs `kinit`).
*Files:* `rsync_push_configs_s2_claude.sh`

### Chunk 5 -- cross-site validation (BLOCKING before bulk launch)
Each remote recomputes **3 configs that LOCAL has already done**, and we compare the h5 against the
local copy. Expected agreement is at the CG tolerance ($10^{-8}$), **not** bit-for-bit -- see section 6.
This is the single check that catches a wrong geometry file, a wrong $\ell$ convention, a wrong
`--spin-dilution` flag, or a mismatched `nu0`.
*Files:* `check_conn_s2_claude.py` (comparison mode)

### Chunk 6 -- bulk run
Each site runs its assigned units (its ensembles x 4 offsets). All resumable and idempotent. Re-run
after each L3/L4 config top-up with the new `--kmax`.
*Files:* `run_conn_s2_claude.sh`, `run_conn_s2_scc_claude.sh`, `conn_stride2_cmds_fnal_claude.txt`

### Chunk 7 -- pull back and merge
Pull `data_*/corr_ylm_conn_t00_nhits1_s1/` from both remotes into the local tree (paths already match,
so it is a straight overlay), then run the completeness check: expected vs present per (ensemble,
offset), `complete` present in every file, no duplicates.
*Files:* `rsync_pull_conn_s2_claude.sh`, `check_conn_s2_claude.py`

### Chunk 8 -- re-analysis
Re-run the conn analyses on the denser grid and re-check the fit windows, which were chosen on
stride-10 data. The L3/L4 windows are currently **provisional copies of L2's** and NM has already
called those fits "terrible", so they must be redone from `effmass_plateau_tables_claude.md` regardless.
*Files:* `effmass_conn_claude.py`, `effmass_plateau_tables_claude.py`,
`effmass_axial_l3_expfit_claude.py`, `axial_ratio_l3_l1_claude.py`

### Chunk 9 -- record
*Files:* `jj_sync_blackboard_claude.md`, memory `project_redo_obs_flow.md`,
`project_jj_fermilab_sync.md`

## 6. Reproducibility across sites

The source is stochastic but the seed is deterministic: `seed_str = ESNID + "_k<k>_h<h>"` hashed to an
`int` by `seed_from_string` (`:379-380`), so every site draws the **same** $Z_2$ wall for a given
$(\text{ensemble}, k, \text{hit})$.

What this does and does not give:
- **Does** give: the same estimator, so results from different sites are statistically interchangeable and
  can be pooled without reweighting.
- **Does not** give: bit-identical h5. Different GPU architecture and different CG iteration counts
  change the last digits; agreement is at `Comp::TOL_OUTER` $= 10^{-8}$. **So validation is by numerical
  comparison, never by md5.**

## 7. Failure modes to guard

| risk | guard |
|---|---|
| two sites compute the same $k$ | disjoint residue classes by construction; `complete`-gating makes an overlap cost a stat, not a solve |
| a site silently has a different geometry (esp. `n3`) | Chunk 5 cross-check against local configs |
| L3/L4 `last` grows mid-run, sites disagree on `--kmax` | `--kmax` passed at launch as `last+1`; re-run tops up |
| partial h5 from a killed job | driver writes `.tmp` then `std::filesystem::rename`, and `complete` is the last dataset written -- a killed job leaves no half-file in place |
| pull-back overwrites good local data | `rsync` without `--delete`; identical paths mean an overlay, and the completeness check runs after |
| ensembles still generating (L3/L4) | fine -- measurement never touches the chain; no fork hazard, unlike the 2026-07-22 HMC handoff |

## 8. Decisions (RESOLVED with NM, 2026-08-10)

**Q1. Pilot first, or commit to stride 2? -> STRIDE 2 NOW.**
Chunk 0 ($\tau_\text{int}$ pilot) is DROPPED. The plan proceeds directly to the full stride-2 grid; the
autocorrelation caveat of section 1 is accepted as a known risk.

**Q2. Split by offset or by ensemble? -> HYBRID: job by offset, machine by ensemble.**
Each ensemble belongs to exactly ONE site; within that site its work is split into the four offset jobs
of section 2.1. Consequences:
- Each site ships only the configs for ITS ensembles, not all 13.6 GB.
- Every ensemble ends up complete at stride 2 regardless of which site is slower, so the $(a/R)^2$
  extrapolation never sees per-ensemble statistics imbalance -- this is precisely what the hybrid buys
  over a pure by-ensemble split.
- The four offsets are the natural intra-site parallelism (4 packed workers LOCAL; 4 queued jobs on SCC).

**Q3. Measurement on the FNAL and SCC allocations? -> YES.**
If an allocation runs out, revisit then.

**Q4. FNAL batch system. -> OUT OF SCOPE FOR LOCAL.**
NM handles FNAL job submission. Deliverable for FNAL is the assignment table plus the exact per-unit
command lines, NOT a batch script. `run_conn_s2_fnal_claude.sh` is REMOVED from section 4.

**Q5. Extend the split to disc? -> NO.**
Conn only. Disc keeps its own separate pending nhits=2 job (`run_disc_h2_claude.sh`).

**Q6. Thermalization cut. -> UNCHANGED.**
$k < 20$ in TRAJECTORY units (so 10 samples on a stride-2 grid, was 2). No analysis change.

## 9. Open questions (still need NM)

**Q7. Relative site capacity.** The by-ensemble machine assignment needs a capacity weight per site.
The assignment generator defaults to **equal thirds** by worker-h; changing that is one line
(`WEIGHTS` at the top of `conn_stride2_assign_claude.py`). LOCAL is 2x TITAN V and is also carrying the
pending disc nhits=2 job. FNAL and SCC GPU counts available for this are unknown here.
