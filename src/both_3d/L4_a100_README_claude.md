# L=4 free current-correlator computation on A100-80 (handoff)

## Goal / context

Part of the `qed3` free-field CFT-verification work (`jj_exact_freefield_impl_plan_claude.md`).  We test
the analytic CFT current-correlator formulas (Eq. 4.31 $G_t$, Eq. 4.28 $G_s$, parameter-free ratio
$G_s/G_t=-(D-1)=-2$) on the **free** (conformal) lattice, refining $L=N_\text{refine}$ to **reduce cutoff
effects**.  The observable here is **method B = the ultralocal ("local") current contracted with the lattice
overlap propagator $P=D_\text{ov}^{-1}$**.  The connected ratio $G_s/G_t$ should approach $-2$ as $L$ grows:

| L | sites | $N=2\,(10L^2{+}2)\,N_t$ | status | $G_s/G_t$ |
|---|---|---|---|---|
| 1 | 12  | 3072  | done (Titan V) | $\approx -12$ (coarse) |
| 2 | 42  | 10752 | done (Titan V) | (in progress) |
| **4** | **162** | **41472** | **THIS JOB (needs A100-80)** | target $\to -2$ |

$N_t=128$.  L=4 dense matrices are $N^2\cdot16 = 27.5$ GB each, so the dense $D_\text{ov}$ build + cuSOLVER LU
need an **A100-80**; a 12 GB Titan V cannot hold them.

## Two stages

1. **Propagator (A100-critical)** — `L4_a100_prop_claude.cu`: builds the dense overlap operator
   $D_\text{ov}$ column-by-column (N `from_cpu` applies), LU-inverts $D_m=D_\text{ov}+m$ (here $m=0$) with
   cuSOLVER `Zgetrf`/`Zgetrs`, and writes `data_free_vmRe0.000000vmIm0.000000/prop_deter_L4/Dinv.0.h5`
   (datasets `Dm_inv` $=D_\text{ov}^{-1}$ and `Dov`; ~55 GB total).  Has a built-in self-check that prints
   `[check] || D_m . D_m^{-1} - I ||_F` (should be ~1e-13).
2. **B = local op + lattice $P$ (host)** — `jj_local_deter_claude.cu` (from the repo): builds the ultralocal
   current $W^P(t)$ (bare $\hat e\cdot\sigma$, no $\Omega$, no $-r\sigma_0$), forms $A(t)=W^P(t)\,P$ and
   $\mathrm{conn}(t_0,t)=\mathrm{tr}(A(t_0)A(t))$ per timeslice, writes
   `.../corr_deter_local_L4/corr.0.h5`.  Host-memory bound (no GPU needed).
3. **D = local op + continuum $G$ (host)** — same `jj_local_deter` binary with
   `--prop-file cont_prop_L4/Dinv.0.h5 --out-tag cont` -> `.../corr_deter_local_cont_L4/corr.0.h5`. Reads the
   continuum propagator instead of the lattice one (same Dinv schema). Host-memory bound; needs
   `cont_prop_L4/Dinv.0.h5` (27.5 GB) present on this server. Stage auto-skips if that file is absent.

`L4_a100_run_claude.sh` runs all three.

## Prerequisites (must exist in the working dir = `src/both_3d/`)

- `includes/` — the full qed3 header set (lattice `s2n_*.h`, `dirac_*.h`, `overlap_wmass_claude.h`,
  `matpoly_claude.h`, `sparse_matrix.h`, `gpu_header.h`, `valence_claude.h`, ...).
- `jj_local_deter_claude.cu` — the stage-2/3 local-op contraction (same dir).
- geometry data `../../geometry/data/*_n4.dat` (pts/links/nns/alpha/omega for $N_\text{refine}=4$).
  **Already generated in the repo** — confirm `alpha_n4.dat`, `omega_n4.dat`, `pts_n4.dat`, ... are present.
- (stage 3 only) `cont_prop_L4/Dinv.0.h5` (27.5 GB continuum propagator) — transfer it to this server if not
  present; stage 3 auto-skips otherwise.
- Eigen 3.4, HighFive, HDF5 1.x.  **Edit `INCLUDES`/`LDFLAGS` in the .sh** to this server's paths.
- Modules: `cuda/12.8`+, `gcc/13.2.0`+.

## Build / arch

- Compile with `-DN_REFINE_CLI=4` (sets $L=4$, $N=41472$; `-DNT_CLI` can override $N_t$, default 128).
- **Change `-arch=sm_70` -> `-arch=sm_80`** (A100) in the .sh `NVCCFLAGS` (done in the provided script).

## Memory

- **GPU (A100-80):** stage 1 holds `d_A` (27.5 GB, $D_\text{ov}$/$D_m$) + `d_B` (27.5 GB, the inverse) +
  getrf workspace + ipiv $\approx 56$ GB.  Fits in 80 GB.  (The column-build phase before LU uses only the
  per-vector overlap buffers, small.)
- **Host RAM:** stage 1 needs ~80 GB (host copy of $D_\text{ov}$ + the inverse real/imag vectors + the Dov
  dump).  Stage 2 (local op) holds $P$ + $A(t_0)$ x n_t0 + $A(t)$ $\approx 4\times27.5 = 110$ GB.
  **Provision ~128-256 GB host RAM.**
  - If host RAM is insufficient for stage 2, the local contraction must be **GPU-ported**: mirror
    `jj_contract_deter_claude.cu::matmul_A` (cuBLAS `Zgemm` for $A=W^P\!\cdot P$) + a GPU trace, streaming
    $A(t)$ over the A100.  The stage-1 propagator is already the hard part and is complete by then.

## Runtime (rough)

L=2 ($N$=10752) dense build was ~98 s + LU ~15 s on a Titan V.  L=4 is $N\!\propto\!3.86\times$ more columns,
each apply heavier, and LU $\propto N^3$ ($\approx 57\times$).  Expect the propagator on the order of
~10-30 min on an A100; stage 2 (host) a few minutes once RAM suffices.  (Times approximate -- please report
actuals.)

## What to verify and return

- Propagator self-check `|| D_m . D_m^{-1} - I ||_F` must be ~1e-13 (printed by stage 1).
- B(L4) connected ratio `Gs/Gt` (printed at the end) should be **closer to $-2$** than L=1 ($\approx-12$)
  and L=2 -- that is the whole point (cutoff reduction).
- **Return:** `corr_deter_local_L4/corr.0.h5` (B) and `corr_deter_local_cont_L4/corr.0.h5` (D), both small,
  under `data_free_vmRe0.000000vmIm0.000000/`.  Keep the 55 GB `prop_deter_L4/Dinv.0.h5` on the server (too
  large to ship; only needed if we later add other contractions).

## Notes / scope

- Only **method B** is requested at L=4.  Method A (exact basis-loop $K$) is L=1-only (too expensive); the
  commutator method was abandoned; method D (local op + a separately-supplied *continuum* propagator) is
  blocked upstream by a continuum-vs-lattice propagator convention mismatch and is NOT part of this job.
- Output schema of `corr.0.h5`: `h0/t0_<b>/{tp,sp}/Vpp/{real,imag}` (length $N_t$, the $1/4\pi$-scaled
  correlator), `h0/disc/{tp,sp}/J` (one-point, ~0), plus `t0s`, `n_t0`, `ls`.
