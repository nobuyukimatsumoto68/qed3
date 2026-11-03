# Incremental distillation-basis extension: L3 Nv 24 -> 32, L4 Nv 48 -> 64 (NM 2026-11-03)

qed3-4e / "Fin: axial sp". NM: "We do this incrementally. L3->32, L4->64" -- add eigenvectors WITHOUT redoing the old
$D_\text{ov}$ solves. Derivation + check on stored data: `final/analysis_axial/peram_extend_nv_note_claude.md`.

**Sources (cite):** distillation, M. Peardon et al., arXiv:0905.2160; the existing driver `distill_peram_mrhs_claude.cu`
(mrhs block CG, fused two-window sources, symmetrized basis BASIS_SYM=1); identity
$(D_\text{ov}^{-1})^\dagger = 1 - D_\text{ov}^{-1}$ from $D_\text{ov}=1+V$, $V$ unitary (paper-21 Eq. A.25, IV.16).

## Physics / goal
Larger basis (finer smearing) gives a second operator size for the variational basis at L3 (sizes 32, 24) and a third
at L4 (64, 48, 24), which fixed L4 l=1 Phi before (48/24). Blocks of $\tau_{ij}(t,t_0)$ with old $=1..N_o$, new
$=N_o{+}1..N_n$:
- old-old: STORED (no work);
- any-new (source new): NEW solves, $N_n-N_o$ sources per source timeslice (8 at L3, 16 at L4), projected on all $N_n$ modes;
- new-old: $\tau_{ij}(t,t_0)=-\tau_{ji}(t_0,t)^*$ ($i\ne j$) from the new solves.
Accuracy of the identity on stored data: ~$10^{-4}$ absolute (not solver precision; the current sp contraction already
relies on it for the backward leg).

## Files
- NEW `src/production/distill_peram_extend_claude.cu`: COPY of `distill_peram_mrhs_claude.cu` (original untouched), changes:
  1. args `--nv-old <No> --old-dir <distill_Nv<No>_sym>`; `--nv <Nn>`; compile `-DNSTACK_CLI=<Nn-No>` (block width = new sources).
  2. basis: same BASIS_SYM dense diagonalization to $N_n$ modes; READ stored `V` (old modes) and USE it (exact match to
     the stored $\tau$); new modes = fresh modes $N_o{+}1..N_n$, re-orthogonalized against the stored ones. Checks logged
     per config: $\max_t\lVert V_\text{stored}^\dagger V_\text{fresh,old}\rVert$ is a unitary (phase/rotation) to 1e-10;
     eigenvalue gap at the cut $\lambda_{N_o+1}-\lambda_{N_o}$ and at $N_n$; stored `evals` reproduced.
  3. solves: only the new sources (fused two-window path as now).
  4. assemble $\tau$ ($N_n\times N_n$) from the three blocks; write `data_<ens>/distill_Nv<Nn>_sym_ext/peram.<k>.h5`
     (NEW dir, never overwrites), same layout + `meta/extended_from` (= old Nv), `meta/identity_block` = 1.
  5. `tau_gw`: OPEN (see below).
- NEW handoff `src/production/run_distill_extend_L3_claude.sh` (local) and `..._L4_claude.sh` (local or SCC): build with
  `/usr/local/cuda-12.9/bin/nvcc -ccbin g++-12`, tee logs, no rm, NM runs.
- analysis side: `interacting_vsh_nvgevp_claude.py` already takes `SIZES`; run with NVDIR = the new dir.

## Chunks
1. Driver copy + edits (items 1-4). Files: `distill_peram_extend_claude.cu`.
2. Validation on ONE config per L (NM runs via `tmp_claude.sh`): (a) extended file vs a FROM-SCRATCH Nv<Nn> run of the
   original driver on the same config: old-old identical, any-new to solver tolerance, new-old to ~1e-4; (b) VSH
   elemental correlators from both files compared (relative difference per dt). Files: `tmp_claude.sh`, logs.
3. Production handoff scripts (L3 490 cfg, L4 390 cfg). Files: `run_distill_extend_L{3,4}_claude.sh`.
4. sp analysis on the new sizes. Files: existing drivers/plotters.

## Cost (from the original runs; linear in the number of new sources)
L3 24->32: 1/3 of the original L3 run (~1.7 days on 1 GPU locally). L4 48->64: 1/3 of the original (~45 L40S-h on SCC,
or ~110 TITAN V GPU-h locally). Disk: L3 ~40 GB, L4 ~120 GB (new full files; old files kept).

## Open questions (NM)
1. `tau_gw` in the extended files: sp uses only `tau`; `tau_gw` users are the FS/two-meson scripts (two-meson dropped
   10-09). Options: omit; or store it for the new-source columns only; or compute the new-old block with one extra
   $D_\text{ov}^\dagger$ apply per old source (no solve).
2. L4: run locally (GPUs busy with the L2 top-up until ~Nov 3-4) or on SCC (`run_distill_sym_L4_scc_claude.sh` pattern)?
3. Validation run before production: OK to spend one from-scratch Nv32 (L3) and Nv64 (L4) config?
