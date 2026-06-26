# L=2 massless local-current $Y_{\ell m}$ tower run -- impl plan

## Goal

Run the bare local-current $Y_{\ell m}$ tower (conn + disc) on the **L=2 massless interacting**
sea ensembles, Nf = 2, 4, 6. Same observable as the L1 massless production
(`run_ylm_prod_claude.sh`), only the lattice refinement (L1 -> L2) and the **config sampling
interval (stride)** change.

Reference thread: memory `project-ylm-tower` (L1 massless production + massive Real-m run).

## What changes vs the L1 production

1. **Lattice size L1 -> L2** = compile-time `N_REFINE` 1 -> 2 (sites $10L^2+2$: 12 -> 42, links 90).
   - The two stoch programs currently hardcode `constexpr int N_REFINE=1;`
     (`jj_local_ylm_conn_stoch_claude.cu:66`, `jj_local_ylm_disc_stoch_claude.cu:64`).
   - Adopt the existing deter-pipeline convention (`jj_kbuild_exact_claude.cu:64`):
     ```
     #ifndef N_REFINE_CLI
     #define N_REFINE_CLI 1
     #endif
       constexpr int N_REFINE=N_REFINE_CLI;
     ```
     Old `constexpr int N_REFINE=1;` line commented out in place (preserve-original style).
   - Compile L2 binaries with `-DN_REFINE_CLI=2`. L1 default unchanged.
   - Geometry `*_n2.dat` confirmed present in `../../geometry/data/`; `Base(2)` loads it.

2. **Ensembles**: `Nf{2,4,6}_gsq8.000000at0.200000nu01.000000nt128L2/` (massless sea).
   - Available checkpoints (k range): Nf2 = 1..1601, Nf4 = 1..518, Nf6 = 1..289.
   - Output auto-tags from `--ens-dir` -> `data_Nf<nf>_...nt128L2_vmRe0.000000vmIm0.000000/`
     `{corr_ylm_conn_t00_nhits1_s1, corr_ylm_disc_tb2_nhits1}/`. No code change for naming.

3. **Config sampling interval (the parameter to change).** L1 used KMIN=40, STRIDE=8, NCONF=140
   (k=40..1152). L2 chains are shorter per Nf, so stride/count are set per ensemble. Thermalization
   skip KMIN=40 kept. Proposed (PENDING user confirmation -- see Open question):

   **DECIDED (user): uniform STRIDE=4, KMIN=40**, NCONF per chain length:

   | Nf | chain | KMIN | STRIDE | NCONF | KMAX(excl) | k range    |
   |----|-------|------|--------|-------|------------|------------|
   | 2  | 1601  | 40   | 4      | 140   | 600        | 40 .. 596  |
   | 4  | 518   | 40   | 4      | 120   | 520        | 40 .. 516  |
   | 6  | 289   | 40   | 4      | 63    | 292        | 40 .. 288  |

   Nf2 capped at 140 (chain has ~390 available at stride 4) to bound run cost; Nf4/Nf6 use the
   full available chain.

4. **Disc dilution `--disc-tblock` = 2 unchanged** (interval $N_t/2=64$; the bias fix is temporal,
   L-independent). DISC_TB stays 2.

## Files

- `jj_local_ylm_conn_stoch_claude.cu` -- add `#ifndef N_REFINE_CLI` guard (minimal, ~3 lines).
- `jj_local_ylm_disc_stoch_claude.cu` -- same guard.
- NEW `run_ylm_prod_L2_claude.sh` -- copy of `run_ylm_prod_claude.sh`: `-DN_REFINE_CLI=2`, L2 ENS(),
  per-Nf KMIN/STRIDE/NCONF, distinct binary names (`*_L2.o`), distinct logs (`ylm_*_L2_nf*_claude.log`).
- NEW `run_ylm_disc_tb2_L2_claude.sh` -- copy of `tmp_ylm_disc_tb2_claude.sh`, same L2 changes.
  (Or fold conn+disc into the single L2 script with both phases, like run_ylm_prod -- TBD.)

## Implementation chunks

- **Chunk 1**: add the `N_REFINE_CLI` guard to both stoch programs. Files: the two `.cu`.
- **Chunk 2**: write the L2 run script(s) with the per-Nf sampling table + `-DN_REFINE_CLI=2`.
  Files: `run_ylm_prod_L2_claude.sh` (+ disc).
- **Chunk 3**: USER runs the script (MPS up); read back logs; verify per-config `.h5` accrue and
  resume-skip works.

## Open question

- **Config sampling stride/count**: confirm the per-Nf table above, or give preferred
  stride / NCONF (decorrelation is the user's call; I can also match L1 exactly where the chain allows).
- GPU/MPS layout: DECIDED (user) -- ALL on GPU 0 (dev=1 not used). Run Nf=2,4 simultaneously on GPU 0
  (2 procs/GPU under MPS, fits the 12 GB TITAN V); Nf=6 run later via
  `bash run_ylm_prod_L2_claude.sh 6`. Script now takes the Nf list as positional args (default "2 4").
