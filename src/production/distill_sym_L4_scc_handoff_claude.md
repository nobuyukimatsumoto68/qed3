# Handoff: symmetrized-basis distillation perambulators for L4 Nf2 gsq2.0 on BU SCC

For a remote agent working on SCC (`scc1.bu.edu`, project `/projectnb/qfe/nmatsum/qed3/src/production`). Author: qed3-ef ("Fin: Two-meson") for NM, 2026-09-22.

## 1. Goal

Produce exact-distillation perambulators for the **L4 Nf2 gsq2.0** massless ensemble using the **symmetrized** distillation basis (`-DBASIS_SYM=1`): $V(t)$ = lowest $N_v$ eigenvectors of the timeslice operator
$$\tilde D^\dagger\tilde D + \tilde D\tilde D^\dagger$$
(one-sided $\tilde D^\dagger\tilde D$ was the old `_v2` basis). Same run parameters as the L1/L2 `_sym` sets already delivered locally.

| item | value |
|---|---|
| ensemble | `Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000` (SCC-generated; HMC may still be appending) |
| lattice | N_REFINE=4, $N_s=162$, $2N_s=324$ (basis is TRUNCATED at $N_v=24$) |
| $N_v$ | 24 (open question 6.1) |
| sources | nsrc=2 fused, tsrc 0 and 64, twin=32, tsrc0=0 |
| stride | 2 over k (open question 6.2) |
| output | `data_<ENS>/distill_Nv24_sym/peram.<k>.h5` |
| source | `distill_peram_mrhs_claude.cu` + `includes/` (unchanged code; only the compile flag differs) |

## 2. Files (all `_claude`, in `src/production`)

Sync these from the local tree into SCC's `src/production` (from the local machine):
```
rsync -av distill_peram_mrhs_claude.cu includes/ \
  tmp_build_distill_sym_L4_scc_claude.sh run_distill_sym_L4_scc_claude.sh run_wrapper_distill_sym_L4_scc_claude.sh \
  nmatsum@scc1.bu.edu:/projectnb/qfe/nmatsum/qed3/src/production/
```
Also needed on SCC: `../../qfe_mod/include` (via `$QED3_INC` from `env.sh`), `../../geometry/data/*_n4.dat`, HighFive at `/projectnb/qfe/nmatsum/opt/highfive/include`, modules `hdf5/1.10.10`, `gsl` (same as the L4 conn jobs -- see `run_wrapper_conn_s2_scc_claude.sh`).

- `tmp_build_distill_sym_L4_scc_claude.sh` -- build-only smoke test, both arches (`sm_70` V100, `sm_80` A100). Flags: `-DN_REFINE_CLI=4 -DNSTACK_CLI=24 -DBASIS_SYM=1`.
- `run_distill_sym_L4_scc_claude.sh` -- SGE batch script = one worker unit (one k-offset on one GPU). Parameters via `qsub -v`.
- `run_wrapper_distill_sym_L4_scc_claude.sh` -- login-node wrapper: build, live `--kmax`, NUNITS staggered units x N_CHAIN `-hold_jid` links.

## 3. Procedure

1. Build smoke test: `bash tmp_build_distill_sym_L4_scc_claude.sh 2>&1 | tee tmp_build_distill_sym_L4_scc_claude.log`. Expect `BUILD OK (both arches)`.
2. Dry run: `DRYRUN=1 bash run_wrapper_distill_sym_L4_scc_claude.sh` -- check the printed `qsub` lines (kmin offsets 1,3,5,7; worker stride 8; kmax = current max ckpoint + 1).
3. Submit: `bash run_wrapper_distill_sym_L4_scc_claude.sh` (defaults: `NV=24 STRIDE=2 NUNITS=4 N_CHAIN=4 H_RT=12:00:00 ARCH=sm_70`). Use `ARCH=sm_80` for A100s if V100 queue is long. Both arches may run concurrently on DIFFERENT units only (same unit twice = same k set; the skip guard prevents overwrite but wastes time).
4. Monitor: `qstat -u nmatsum`; `ls data_<ENS>/distill_Nv24_sym/peram.*.h5 | wc -l`; per-job logs `dsym_L4g2.0_u<u>_c<c>.o<jid>`.
5. Re-submit is always safe: the binary skips any `peram.<k>.h5` that exists (`distill_peram_mrhs_claude.cu:471-478`).

## 4. Validation (when done, or on a partial set)

Run on SCC (python3 + h5py + numpy):
```python
import glob
import numpy as np
import h5py
d = "data_Nf2_gsq2.000000at0.200000nu01.000000mRe0.000000mIm0.000000nt128L4_hb0.400000-1.000000/distill_Nv24_sym/"
fs = sorted(glob.glob(d + "peram.*.h5"), key=lambda s: int(s.split(".")[-2]))
print("n perams", len(fs))
for fn in fs[:2] + fs[-1:]:
    h = h5py.File(fn, "r")
    tau = np.array(h["peram/tau/real"]) + 1j * np.array(h["peram/tau/imag"])
    ev = np.array(h["evals"])
    h.close()
    dev = max(np.abs(np.diag(tau[s, t, t]).real - 0.5).max() for s in range(tau.shape[0]) for t in range(tau.shape[1]))
    print(fn, tau.shape, "max|Re diag tau - 1/2| =", dev, "finite", np.isfinite(tau).all(), "evals ascending", bool(np.all(np.diff(ev, axis=1) >= -1e-12)))
```
Pass criteria (matching the delivered L1/L2 sets): shape `(2,32,32,24,24)`; `max|Re diag tau - 1/2|` about 1e-5 (= `overlap_tol`, the GW identity $\mathrm{Re}\,\tau_{ii}=1/2$); all finite; evals ascending; uniform file sizes (`stat -c %s ... | sort | uniq -c`); no gaps on the stride grid.

## 5. Sync back to local

From the local machine, extend the pattern of `rsync_scc_claude.sh`:
```
rsync -avz --prune-empty-dirs --include='data_Nf2_gsq2.0*L4_hb*/' --include='distill_Nv24_sym/' --include='peram.*.h5' --exclude='*' \
  nmatsum@scc1.bu.edu:/projectnb/qfe/nmatsum/qed3/src/production/ /mnt/barracuda22/qed3/qed3/src/production/
```
Local consumers select the set with `NVDIR=distill_Nv24_sym`, `LREF=4`; use `tau` only (not `tau_gw`) and `MODE_CONTACT=1` (mode-space contact) for any truncated-basis contraction.

## 6. Open questions for NM

1. **$N_v$**: 24 keeps the same basis size as L1/L2 (24 of 324 at L4 = heavily truncated). Larger $N_v$ (e.g. 48) costs ~linearly in solves and h5 size; `NV=<n>` on the wrapper handles it (binary is rebuilt with `NSTACK_CLI=<n>`).
2. **Stride**: 2 gives ~400 perams from ~800 configs; stride 1 doubles cost.
3. **Arch**: V100 (sm_70) default; A100 (sm_80) also strong-FP64 and fine.

## 7. Cost estimate

L2 `_sym` (2Ns=84) ran at 8.2 min/peram/worker with 2 MPS workers per TITAN V. Scaling by sites (162/42 ~ 3.9x) and assuming one job per GPU: ~15-25 min/peram on a V100. ~400 perams / 4 units ~ 100 per unit ~ 25-40 h per unit -> 3-4 chained 12 h links. Wall time ~1.5-2 days if 4 GPUs are available.

## 8. Rules

No `rm`, no killing jobs, never write into `distill_Nv24_v2` / production dirs (`--out-suffix _sym` is fixed in the scripts). Measurement only: reads `ckpoint_lat.*`, never modifies configs.

## 9. Refs

Distillation: Peardon et al. arXiv:0905.2160. Basis/regen background: `peram_regen_changelog_claude.md`; local L1/L2 runner `run_distill_sym_L1L2_Nf2_claude.sh`.
