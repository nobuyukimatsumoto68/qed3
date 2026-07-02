# Handoff: running the JJ (ylm conn+disc) calculation on REMOTE (Fermilab)

For the REMOTE agent on `nmatsum@lq.fnal.gov`. End-to-end: **compilation -> job submission**.
Companion docs in this dir: `jj_sync_blackboard_claude.md` (comm log, Q1-Q5), `jj_ensembles_claude.txt`
(manifest), `fermilab_dirs.txt` (canonical config-dir list), `fermilab_ncfg_claude.txt` (config counts).

Local reference drivers you are mirroring: `run_ylm_prod_L2_massive_claude.sh` (recommended base,
auto-detects config sampling), `run_ylm_prod_claude.sh` (L1), `run_ylm_massive_pair_claude.sh` (L1
massive, MPS pack).

---

## 0. Scope (Q1, resolved)

Compute the local-current $Y_{\ell m}$ tower (**conn + disc**) on the jj-first ensembles:

- **Priority = gsq8 L1 + L2**, including the `heavy/` L1/L2 masses (all COMPLETE @319 cfg).
- The genuine **compute gaps** (no jj exists anywhere -- LOCAL or REMOTE): **`heavy` L1/L2 @319**
  (ready NOW), all **L4** (Family-B @119, heavy L4 ~58-66, still growing), **gsq12** (1-5 cfg, tiny).
- **DEFER** L4 and gsq12 until they accumulate.
- **OUT of scope:** the parity `mIm!=0` set -- the stochastic programs `assert(!parity)` and no
  analysis is planned on it. Only `mIm=0` (massless or real $m_F$) here.

Valence mass = **sea mass** (unitary point) = the `mRe` in the ensemble dir name; `mIm=0`.

---

## 1. The two compute programs

- `jj_local_ylm_conn_stoch_claude.cu` -- gauge-connected per-$m$ vector+axial tower, single $t_0=0$
  wall source, spin dilution.
- `jj_local_ylm_disc_stoch_claude.cu` -- disconnected one-point $J^a_{\ell m}$ + $\sigma_\text{PS}$,
  time+spin dilution, `--disc-tblock 2` (bias fix; interval $N_t/2=64$).

Physical tower (LOCAL analysis) = $-C_\text{conn}+C_\text{disc}$. You only PRODUCE the two h5 trees;
analysis stays on LOCAL.

---

## 2. Prerequisites (check before building)

1. **Code**: this repo on branch **`develop`** at the HEAD NM pins in the blackboard (Q5). The `.cu`
   above + `includes/` (esp. `valence_claude.h`, `s2n_simp.h`, `sparse_dirac_claude.h`,
   `matpoly_claude.h`, `gauge_ext.h`, `action_ext.h`).
2. **Working directory MUST be `<repo>/src/both_3d/`.** The geometry path is HARDCODED as
   `"../../geometry/data/"` (`jj_local_ylm_conn_stoch_claude.cu:84`), and output `data_*` is written
   relative to CWD. So `cd .../src/both_3d` before running; `../../geometry/data/` must resolve.
3. **Geometry data** for the refinement you build: `pts_n<L>.dat`, `pts_dual_n<L>.dat`,
   `nns_dual_n<L>.dat`, `links_n<L>.dat`, `dual_links_n<L>.dat`, `face_n<L>.dat` under
   `../../geometry/data/` (L1 -> `_n1`, L2 -> `_n2`). Confirm present.
4. **Configs synced** on `/lustre2` at the abs paths in `fermilab_dirs.txt` (you have these).
5. **Deps**: CUDA for **A100 = `sm_80`** (NOT `sm_70`), Eigen 3.4, HighFive + HDF5, GSL, OpenMP.
   `source /home/nmatsum/env.sh` (as `run_nf_fermilab_affine_claude.sh` does) should set module env;
   adjust the `-I`/`-L` include/lib paths below to the cluster's actual locations.
6. **Optional MPS** for packing several small jobs on one A100:
   `nvidia-cuda-mps-control -d` (each proc is only ~400 MiB, so an A100 packs many).

---

## 3. Compilation

**One binary PER refinement L** (N_REFINE is compile-time; you cannot switch L at runtime). Build
`_L1.o` and `_L2.o` pairs. Change `sm_70 -> sm_80` for A100 and fix include/lib paths for the cluster:

```
NVCC=nvcc
# A100 -> sm_80 ; L via -DN_REFINE_CLI=<L>
NVCCFLAGS="-arch=sm_80 -g -O3 -std=c++20 -DN_REFINE_CLI=<L> -lcublas -lcusolver -lcusparse -lgomp -Xcompiler -fopenmp"
# adjust these to the cluster (env.sh likely provides EIGEN/HIGHFIVE/HDF5 roots):
INCLUDES="-I./includes/ -I<EIGEN_ROOT> -I<HIGHFIVE_ROOT>/include -I<HDF5_ROOT>/include"
LDFLAGS="-L<HDF5_ROOT>/lib -lhdf5 -lgsl -lgslcblas -lm"

# L1 (-DN_REFINE_CLI=1):
$NVCC ${NVCCFLAGS/<L>/1} $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o jj_local_ylm_conn_stoch_L1.o
$NVCC ${NVCCFLAGS/<L>/1} $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o jj_local_ylm_disc_stoch_L1.o
# L2 (-DN_REFINE_CLI=2):
$NVCC ${NVCCFLAGS/<L>/2} $INCLUDES $LDFLAGS jj_local_ylm_conn_stoch_claude.cu -o jj_local_ylm_conn_stoch_L2.o
$NVCC ${NVCCFLAGS/<L>/2} $INCLUDES $LDFLAGS jj_local_ylm_disc_stoch_claude.cu -o jj_local_ylm_disc_stoch_L2.o
```

(LOCAL uses `-arch=sm_70` on TITAN V; the rest of the flag string is identical to
`run_ylm_prod_L2_massive_claude.sh`.)

---

## 4. Program invocation (CLI)

Flags (parsed in the `.cu`): `--ens-dir <path>/` `--kmin` `--stride` `--kmax` `--nhits`
`--t0 <n>` `--spin-dilution` (conn) `--disc-tblock <tb>` `--time-dilution <n>` (disc)
`--mass-re <m>` `--mass-im <m>` (omit both for massless).

- **conn**: `--t0 0 --spin-dilution`
- **disc**: `--disc-tblock 2`
- **`--ens-dir` = the ABSOLUTE `/lustre2/...` path** from `fermilab_dirs.txt` (with trailing `/`).
  The program reads `<ens-dir>/ckpoint_lat.<k>` (`...conn...:293`). LOCAL uses a basename because
  configs are local; on REMOTE pass the full path.
- **`--mass-re` = EXACT full-precision sea mass** (unitary), NOT the `%f`-rounded dir tag. E.g. the
  Family-B L2 masses are `0.0105724914705, 0.0528624573524, 0.1057249147049, 0.2114498294097`. For
  the `heavy` set (`mRe0.422900/0.845799/1.268699`) pass the exact HMC sea masses you generated with.
  Massless -> omit `--mass-re` entirely. **`--mass-im 0` always (no parity).**
- **Output** (Q4, the `.cu`'s native naming, do NOT invent a tree): written to CWD as
  `data_<ESNID>/corr_ylm_conn_t0<t0>_nhits<H>_s<0|1>/corr.<k>.h<h>.h5` and
  `data_<ESNID>/corr_ylm_disc_tb<tb>_nhits<H>/corr.<k>.h<h>.h5`, where
  **`ESNID = <ens_basename>_vmRe<mass_re>vmIm<mass_im>`** (`...:271`). Per-config atomic `.h5` +
  resume-skip (safe to re-run; completed configs are skipped).

Example (one L2 massive ensemble, conn then disc):

```
cd <repo>/src/both_3d
ENS=/lustre2/affine/Nf4_gsq8.000000at0.200000nu01.000000mRe0.105725mIm0.000000nt128L2
M=0.1057249147049
CUDA_VISIBLE_DEVICES=0 ./jj_local_ylm_conn_stoch_L2.o --ens-dir "$ENS/" --kmin 40 --stride 4 --kmax 320 --nhits 1 --t0 0 --spin-dilution --mass-re $M --mass-im 0
CUDA_VISIBLE_DEVICES=0 ./jj_local_ylm_disc_stoch_L2.o --ens-dir "$ENS/" --kmin 40 --stride 4 --kmax 320 --nhits 1 --disc-tblock 2 --mass-re $M --mass-im 0
```

---

## 5. Config sampling (Q2: provisional; choose at jj time)

Convention (mirror the reference driver): **`kmin = first k >= 40`** (thermalization skip),
**`stride = SAMP * base`** where `base` = spacing of the first two `ckpoint_lat.<k>` (contiguous
storage -> base=1), **`kmax = last + 1`**. LOCAL used `SAMP=8` at L1 (stride 8) and `SAMP=4` at L2/L4
(stride 4). These caps are NOT fixed -- pick the sample count you want; `src_ncfg` in the manifest is
the durable number. The programs **BREAK at the first missing config**, so the grid must hit existing
files -- auto-detect from the actual `ls ckpoint_lat.*` (see `run_job` in the reference driver).

---

## 6. Job submission (SLURM)

Base the sbatch header on `run_nf_fermilab_affine_claude.sh`:

```
#!/bin/sh
#SBATCH --account=affine.lq2_gpu
#SBATCH --qos=normal
#SBATCH --partition=lq2_gpu
#SBATCH --gpus=a100:1
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm_%x_%j.out
source /home/nmatsum/env.sh
cd <repo>/src/both_3d
export OMP_NUM_THREADS=4
# ... run conn + disc for the assigned ensemble(s) ...
```

Two viable models:
- **(a) one sbatch job per ensemble** (runs conn then disc serially). Simplest; scales by array/loop
  over the jj-first ensembles.
- **(b) MPS-packed** (like the local 2-wide pool in `run_ylm_prod_L2_massive_claude.sh`): start
  `nvidia-cuda-mps-control -d`, run several conn/disc procs concurrently on the one A100 (memory is
  not the limit; latency/compute-bound). Good for the many small L2 ensembles.

**Driver skeleton** (adapt paths; reads `fermilab_dirs.txt`, filters to jj-first, auto-detects grid):

```
#!/usr/bin/env bash
set -u
cd <repo>/src/both_3d
KMIN_TRAJ=40
SAMP_L1=8
SAMP_L2=4
NHITS=1
# jj-first filter: gsq8 L1/L2 incl heavy L1/L2; skip L4, gsq12, mIm!=0
grep -vE '^\s*#|^\s*$' fermilab_dirs.txt \
  | grep -E 'gsq8\..*nt128L[12]$' \
  | grep -vE 'mIm0\.(01|05|1|2)0000[0-9]' \
  | while read -r ENS; do
      base=$(basename "$ENS")
      L=$(printf '%s' "$base" | sed -E 's/.*nt[0-9]+L([0-9]+)$/\1/')
      case "$L" in 1) BIN=_L1.o; SAMP=$SAMP_L1;; 2) BIN=_L2.o; SAMP=$SAMP_L2;; esac
      # exact sea mass: look up the full-precision value for this dir (dir tag is %f-rounded);
      # massless (no mRe in name) -> MASS empty
      MRE=$(printf '%s' "$base" | grep -oE 'mRe[0-9.]+' | sed 's/mRe//')
      MASS=""   # <-- REMOTE: map $MRE (rounded) -> exact sea mass here; leave empty if massless
      # auto-detect kmin/stride/kmax from the on-disk configs:
      mapfile -t ks < <(ls "$ENS"/ckpoint_lat.* 2>/dev/null | sed 's#.*ckpoint_lat\.##' | grep -E '^[0-9]+$' | sort -n)
      [ "${#ks[@]}" -eq 0 ] && { echo "SKIP $base (no configs)"; continue; }
      base_sp=$(( ${#ks[@]} >= 2 ? ks[1]-ks[0] : 1 )); [ "$base_sp" -le 0 ] && base_sp=1
      kmin=""; for k in "${ks[@]}"; do [ "$k" -ge "$KMIN_TRAJ" ] && { kmin=$k; break; }; done
      [ -z "$kmin" ] && kmin=${ks[0]}
      stride=$(( SAMP * base_sp )); kmax=$(( ks[-1] + 1 ))
      MASSF=(); [ -n "$MASS" ] && MASSF=(--mass-re "$MASS" --mass-im 0)
      # submit (model a: one job per ensemble, conn then disc):
      sbatch --job-name="jj_${base}" --wrap="source /home/nmatsum/env.sh; cd <repo>/src/both_3d; \
        ./jj_local_ylm_conn_stoch${BIN} --ens-dir '$ENS/' --kmin $kmin --stride $stride --kmax $kmax --nhits $NHITS --t0 0 --spin-dilution ${MASSF[*]}; \
        ./jj_local_ylm_disc_stoch${BIN} --ens-dir '$ENS/' --kmin $kmin --stride $stride --kmax $kmax --nhits $NHITS --disc-tblock 2 ${MASSF[*]}" \
        --account=affine.lq2_gpu --qos=normal --partition=lq2_gpu --gpus=a100:1 --time=6:00:00 --cpus-per-task=16 --output=slurm_%x_%j.out
    done
```

The `MASS=""` line is the **one thing REMOTE must fill**: the exact full-precision sea mass per
`mRe` tag (you have these from generation). Everything else auto-derives.

---

## 7. Sanity check FIRST (before fanning out)

Run ONE config on one ensemble (`--kmin K --stride 1 --kmax K+1`), confirm
`data_<ESNID>/corr_ylm_conn_.../corr.<K>.h0.h5` appears and opens (expected groups: the per-$\ell m$
vector+axial tower; disc file has `sigma_PS`/condensate group). Then submit the full set.

---

## 8. Output location + pull-back

Output `data_*` lands in `<repo>/src/both_3d/` (CWD). LOCAL will pull it with a forthcoming
**jj pull-back rsync** that globs `data_*/corr_ylm_{conn,disc}_*/corr.*.h5` from this dir. Keep the
run CWD stable so that glob is trivial. (If you must run from `/lustre2` scratch instead, tell LOCAL
the base dir on the blackboard so the pull-back can target it.)

---

## 9. Report back (blackboard)

When jobs land, append a **Log** entry to `jj_sync_blackboard_claude.md`: which ensembles produced
conn+disc, config counts, the **commit built**, and the **output base dir**. Flag any build/runtime
issue (missing geometry `_n<L>.dat`, HDF5/Eigen path, config gap that tripped the BREAK) so LOCAL can
adjust. NM iterates the pinned HEAD as needed (Q5).
