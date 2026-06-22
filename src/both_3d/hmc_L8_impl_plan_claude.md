# L=8 refinement: implementation plan

Add an $L=8$ (N_REFINE=8) stream to the diagonal-mass QED3 HMC campaign, run under the
**affine** account alongside the existing $L=2$. Diagonal measure-weighted mass
(arXiv:2510.03085); same physical masses at every $L$. Integrator MN2; Zolotarev sign
function for the overlap. See [[l2l4-refinement-run]], [[mass-measure-factor]].

## Goal / physics summary
Pure spatial refinement of the simplicial $S^2\times R$ lattice: the spatial spacing
$a_\text{lat}=$ `mean_ell` halves each refinement ($L1\,1.1071 \to L2\,0.5909 \to L4\,0.2995
\to L8\approx 0.150$), the physical sphere radius $R=1$ fixed. Physical mass is the SAME at
every $L$ (the diagonal factor $A_y/\bar a_s$ carries the $L$-dependence), so masses and
pairs are UNCHANGED from $L2/L4$:
$$m \in \{0.0105724914705,\ 0.0528624573524,\ 0.1057249147049,\ 0.2114498294097\}$$
MPS-packed, two pairs: A = (light 0.010572, heavy 0.211450), B = (0.052862, 0.105725).
`mass_coeff = m * mean_dual_area/mean_ell` depends on GEOMETRY only (NOT on `at`), so it is
fixed by $N\_REFINE=8$ regardless of the `at` choice below.

## BLOCKER: admissibility forces a smaller `at` at L8 (needs a decision)
The driver asserts (overlap kernel admissibility)
$$\sqrt{3}\,\frac{\text{mean\_ell}}{at} - \frac{4}{\sqrt 3} > 0
\quad\Longleftrightarrow\quad \frac{\text{mean\_ell}}{at} > \frac43 .$$
- L2: $0.5909/0.2 = 2.96$ OK; L4: $0.2995/0.2 = 1.50$ OK (tight);
- **L8 @ at=0.2: $0.150/0.2 = 0.75 < 1.333$ -> assert FAILS, binary aborts.**

Max admissible `at` at L8:
$$at_\text{max} = \tfrac34\,\text{mean\_ell}(L8).$$
$\text{mean\_ell}(L8)$ is computed by the binary (printed `# alat`); estimated from the
measured halving trend (L1 1.10715, L2 0.59095, L4 0.29947; ratios 0.534, 0.507 $\to$ 0.5)
and the flat-equilateral side $s\approx 1.2046/N\_REFINE$ (matches measured to $<1\%$ by L4)
as $\text{mean\_ell}(L8)\approx 1.2046/8\times0.997 \approx 0.1503$. Hence
$$at_\text{max} \approx 0.75\times0.1503 \approx 0.1127 .$$
**Decided: at = 0.1** ($\text{mean\_ell}/at = 1.50$, same margin as L4; ample headroom vs the
0.1127 bound even if the true alat is ~0.148). nsteps=16, npole=13 (confirmed by NM).
The natural choice that mirrors the spatial halving is **at = 0.1** (then
$\text{mean\_ell}/at = 1.50$, the SAME margin as L4). But `at` was held FIXED at 0.2 across
L1/L2/L4, so shrinking it changes the temporal discretization. Two sub-choices for `Nt`:

| Option | at | Nt | Nt*at (temporal box) | admiss. | cost vs L4 | notes |
|---|---|---|---|---|---|---|
| **A** keep box | 0.1 | 256 | 25.6 (= L1/L2/L4) | 1.50 | ~8x (4x spatial x 2x Nt) | clean LCP, temporal box preserved; very expensive |
| **B** keep Nt | 0.1 | 128 | 12.8 (half) | 1.50 | ~4x (spatial only) | cheaper, but temporal box differs from L2/L4 -> not a clean temporal LCP |
| C keep at | 0.2 | 128 | 25.6 | 0.75 FAIL | -- | NOT viable (assert aborts; relaxing it breaks overlap gap) |

**RESOLVED (NM, 2026-06-22): at = 0.1, Nt = 128** (Option B). Nt NOT doubled -> temporal box
$Nt\cdot at = 12.8$ is HALF of L1/L2/L4's 25.6 (spatial-continuum target only, not a temporal
LCP). Cost ~8x L4/traj (4x spatial x 2x nsteps). n = 13, nsteps = 16 confirmed (below).

## Extrapolation of n (Zolotarev poles) and nsteps
Both depend on the resolved `at`. Numbers below assume **at = 0.1** (Option A or B).

### nsteps (tmax = 1.0 fixed; MN2)
$nsteps$ history at tmax=1.0: L2 = 5 ($dt=0.20$), L4 = 8 ($dt=0.125$; bumped 6->7->8 to tame
a worst-conditioned-mass blow-up). The fermion force stiffens as the spacing shrinks; the
stable $dt$ scales roughly with the lattice spacing. L2->L4 used $dt$ ratio $0.125/0.20=0.63$
(nsteps ratio 1.6) for a spacing ratio ~0.5. Applying the same factor L4->L8:
$$nsteps_{L8} \approx 8 \times 1.6 \approx 13;\quad \text{2x-(spacing) bound} \to 16.$$
Given L4 already needed bumping for stability and L8 recovery is very expensive, **propose
$nsteps_{L8} = 16$** (conservative; $dt=0.0625$). Adjustable down if first-traj $|dH|$ is small.

### n = npole (Zolotarev), must be ODD
$npole$ history (at=0.2): L2 = 17 (kernel ratio $\lambda_\min/\lambda_\max=0.149$), L4 = 13
(0.210) -- fewer poles as the kernel conditioned better with fixed `at`. BUT shrinking `at` to
0.1 changes the kernel spectrum (worsens conditioning), so the downward L2->L4 trend does NOT
simply continue. Safer to NOT lower below L4. **Propose $n_{L8} = 13$** as the starting guess
and CONFIRM from the printed `# delta` and `# min max ratio` at startup (binary prints both);
bump to 15 if $\delta \gtrsim 10^{-4}$ or the ratio falls outside the Zolotarev window. n must
stay ODD.

**Both n and nsteps are first-run guesses to be confirmed from the startup banner / first-traj
$|dH|$, exactly as L2/L4 were tuned.**

## Files to create / modify
1. **`hmc_fermilab_wmass_L8_claude.cu`** -- copy of `_L4_claude.cu`; change
   `N_REFINE=4 -> 8` (`:56`), `at 0.2 -> 0.1` (`:176`), `Nt` stays 128 (`:60`),
   `Fermion D(DW,mass,13)` npole (`:187`), `nsteps -> 16` (`:259-262`), keep `tmax=1.0`,
   `kmax=80` (`:218`), keep `#define GRAD_L4`. N_SITES=$10\cdot8^2+2=642$, N_LINKS=$30\cdot8^2
   =1920$ follow automatically.
2. **Build**: add L8 target to `rebuild_prod_diag_claude.sh` (generic `%.o:%.cu`,
   `-arch=sm_80`); produces `hmc_fermilab_wmass_L8_claude.o`.
3. **`/lustre2/affine/run_nf_fermilab_mps_L8_claude.sh`** -- copy of the L2 run script;
   `L=8`, `BIN=.../hmc_fermilab_wmass_L8_claude.o`, `--time=16:00:00`, `WALL_SEC=57600`
   (MAX_SEC/BACKSTOP derive). Same MPS daemon + 2 backgrounded clients + 6th-arg MAX_SEC.
4. **`/lustre2/affine/run_wrapper_nf_fermilab_mps_L8_claude.sh`** -- copy of the L2 wrapper;
   `L=8`, points at the L8 run script, SAME masses/pairs, `for Nf in 2 4 6`, afterany chain
   keyed on `hmc_Nf<Nf>_L8_pair<A|B>` via `sbatch --parsable`. dir3 encodes L8 -> no collision
   with L2 dirs in the shared /lustre2/affine cwd.

## STATUS: all files DONE (2026-06-22) -- UNPACKED, split 2+2 across accounts
Driver `hmc_fermilab_wmass_L8_claude.cu`: N_REFINE=8, at=0.1, **Nt=128**, npole=13, nsteps=16,
tmax=1.0, kmax=80, GRAD_L4 (the optimized force codepath, NOT L-specific -- L8 uses it too),
lat every conf (`k_ckpoint=1`), **rng every 2 confs (`k_ckpoint_rng=2`)**, NSTREAMS=1 (unchanged).
Build target added to `rebuild_prod_diag_claude.sh`.

### GPU: UNPACKED (no MPS), one mass per A100
MPS packing benefit declines with L (L2 1.87x -> L4 1.45x) because a bigger local problem
saturates the A100; L8 is 4x L4 (N_LINKS 1920, Nt=128) so a single client fills the GPU ->
packing no longer pays. Run one mass/GPU. (NSTREAMS kept at 1 per NM -- not raised to 4.)

### Four masses SPLIT 2+2 across accounts (Nf=2)
- **qed3** runs the two LIGHTEST: 0.0105724914705, 0.0528624573524.
- **affine** runs the two HEAVIEST: 0.1057249147049, 0.2114498294097.
Spreads compute (both allocations) and disk. NON-MPS single-mass scripts (each account has its
own copy): `run_nf_fermilab_L8_claude.sh` + `run_wrapper_nf_fermilab_L8_claude.sh` in
`/lustre2/affine` and `/lustre2/qed3` (account + MASSES differ); all 4 `bash -n` clean, chmod +x.
Job name `hmc_Nf2_L8_mR<mass>`; afterany chain via `sbatch --parsable`.
OBSOLETE (superseded, hand-off rm to NM): `/lustre2/affine/run_{,wrapper_}nf_fermilab_mps_L8_claude.sh`.

### Disk: rng every 2 confs, split 2+2
Measured rng: L4 0.533 GB, L2 0.134 GB (rng $\propto N\_SITES\cdot Nt$). L8 @ Nt=128
~ $0.533\times3.96 = 2.1$ GB/file. rng every 2 confs -> ~40 kept/stream; 2 streams/account x
kmax=80 ~ **169 GB/account**, well under affine 1 TB / qed3 2 TB (lat every conf but tiny, ~5 MB).

### LAUNCH SEQUENCE
1. `cd /project/qed3/qed3/src/both_3d && bash rebuild_prod_diag_claude.sh` -- builds L2/L4/L8
   (must end `ALL OK`; produces `hmc_fermilab_wmass_L8_claude.o`, shared by both accounts).
2. `cd /lustre2/qed3   && bash run_wrapper_nf_fermilab_L8_claude.sh 1`   (two lightest, qed3).
3. `cd /lustre2/affine && bash run_wrapper_nf_fermilab_L8_claude.sh 1`   (two heaviest, affine).
   Start NCHAIN=1 each; watch the first run before queueing a chain.

### FIRST-RUN CHECKS (critical at a brand-new L / at)
- startup `# alat` = mean_ell(L8) (~0.150 expected) and the admissibility assert PASSES
  (mean_ell/at ~ 1.5; binary would abort otherwise).
- `# delta` and `# min max ratio` -> confirm npole=13 gives delta <~ 1e-4; bump to 15 (ODD) if not.
- first-traj |dH| / acceptance at nsteps=16; sec/traj vs the 16h budget (~4-5 h/traj -> ~3 traj/slot)
  -- if one traj risks exceeding 16h it gets killed mid-traj with NO checkpoint -> chain stalls; then
  reduce nsteps or split differently.

## Remaining open items
1. kmax=80 to start (user-set). Extend later if needed.
2. Cost: L8 draws from the affine allocation SHARED with running L2 -- check `lquota` before a
   long chain.

## Cost / allocation caveat
L8 (Nt=128) draws from the affine allocation SHARED with the running L2 streams. ~8x L4/traj
(4x spatial x 2x nsteps) plus the 16h slot means only ~3 configs/job. Flag the
affine `lquota` (compute AND disk) before committing a long chain (binding constraint per
[[l2l4-refinement-run]]).
