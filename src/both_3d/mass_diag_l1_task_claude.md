# TASK: diagonal measure-weighted mass — implement + L=1 check (local run)

Self-contained brief for the local environment. Full rationale/physics:
`mass_measure_factor_impl_plan_claude.md` (same dir). Source: arXiv:2510.03085.

## STATUS (2026-06-19)
- **Step 1 DONE**: `includes/overlap_wmass_obsolete_claude.h` = frozen scalar reference,
  wrapped in `namespace obsolete`.
- **Step 2 DONE (non-ms path)**: `includes/overlap_wmass_claude.h` now has complex
  diagonal `m_L`: members `M_mass`/`mass_coeff`(Complex)/`d_mtmp`; ctor builds them;
  helpers `apply_mL` / `apply_mLdag`; `mult`/`adj`/`DHD_deviceAsyncLaunch` converted.
- **TODO (do NOT use in production until done)**: the `_ms` variants (`mult_ms` ~:474,
  `adj_ms` ~:497, `DHD_..._ms` ~:559-563) and the **mass inner products** still use the
  scalar `mass` — convert to `apply_mL`/`apply_mLdag` before any HMC run. The L=1 check
  below only exercises the (converted) non-ms `mult`/`adj`/`DHD`.
- **Step 3 WRITTEN**: `test_diag_mass_l1_claude.cu` (simplicial; loops m = 0, 0.1, 0.1 i).
  Parametric in L: build `-DLREF=1|2|4`. **L=1**: operator obsolete-vs-production on BOTH
  `mult/adj/DHD` and `_ms` (~1e-13) + force obsolete-vs-production over all links (~1e-9) +
  force-vs-FD. **L>1**: force-vs-FD only (no scalar reference). Both grad routines via build:
  default (converted) vs `-DGRAD_L4` (block; convert grad_l4 first). Force check mimics
  `hmc_fermilab_claude.cu:357-487`.
- **Remaining for local env**: add a Makefile.fnal target + build + run on GPU (Step 4).
  Check first: the geometry path (`/project/affine/.../geometry/`) and that `Gauge U(base)`
  is a valid cold gauge in your env; bump `Comp::Nt` if Nt=2 is unhappy. Adjust if the
  obsolete-header `namespace obsolete` symbols or include order need tweaks.

## Goal

Replace the scalar fermion mass $D_\text{ov}+m\,\mathbb{1}$ with the measure-weighted
diagonal $D_m=D_\text{ov}+m_L$, $m_L=\text{diag}(m\,A_y/\bar a_s)$ ($m$ = physical mass,
$R=1$ units; $A_y$ = dual-cell area; $\bar a_s$ = `mean_ell`). Then verify that at $L=1$
(uniform measure) it reduces **exactly** to the old scalar operator on a random vector.

## Key facts (already pinned)

- $m_L = \texttt{mass\_coeff}\times\texttt{volume\_matrix}(1)$, where
  `volume_matrix(.,1)` = $\text{diag}(A_y/\bar A)$ (`dirac_ext.h:350`) and
  `mass_coeff = m_phys * lattice.mean_dual_area / lattice.mean_ell`. At $L=1$,
  `volume_matrix(1)`=identity $\Rightarrow m_L=\texttt{mass\_coeff}\cdot\mathbb{1}$.
- Geometry reachable as `DW.lattice.{dual_areas[ix], mean_dual_area, mean_ell}`
  (`s2n_simp.h:18,38,39`).
- Squared operator (real diagonal $M$, uses GW $D^\dagger D=D+D^\dagger$):
  $D_m^\dagger D_m=(1+M)D+D^\dagger(1+M)+M^2$, computed as
  $(1+M)(D_\text{ov}v)+D_\text{ov}^\dagger[(1+M)v]+M^2v$ — same overlap-apply count as the
  old scalar identity (NO extra CG cost).
- GEVP code to mimic for the COO+MatPoly diagonal apply: `saved_scripts/eig.cu:263-286`.
- COO/CSR/MatPoly: `includes/matpoly.h` (`push_back`:81, `on_gpu`:106), `sparse_matrix.h`.

## Steps

### 1. Freeze the scalar reference
```
cp includes/overlap_wmass_claude.h includes/overlap_wmass_obsolete_claude.h
```
In the obsolete copy, wrap all file-scope structs in a namespace so a single test TU can
include both headers without clashes: after the `#include`s (after line 4) add
`namespace obsolete {`, and add a matching `}` at EOF. (Covers `Zolotarev`,
the rational-function struct ~line 117, and `OverlapWMass` → `obsolete::OverlapWMass`.)
Do NOT otherwise edit the obsolete copy — it is the scalar ground truth.

### 2. Edit the production `includes/overlap_wmass_claude.h` (diagonal mass)
- **ctor (`:173-192`)**: keep the `mass` arg but treat it as the physical $m$ (real for our
  runs). Build a member `COO<N> M_mass; DW.volume_matrix(M_mass.en, 1.0); M_mass.do_it();`
  and store `double mass_coeff = m_phys * DW.lattice.mean_dual_area / DW.lattice.mean_ell;`.
  Allocate one scratch device buffer `CuC* d_mtmp` (size `N*CD`).
- **`mult` (`:373`)**: replace `Taxpy_gen(d_res, cplx(mass), d_xi, d_res)` with
  `+ m_L v`: `MatPoly Op_m; Op_m.push_back(cplx(mass_coeff), {&M_mass}); Op_m.on_gpu<N>(d_mtmp, d_xi);`
  then `Taxpy_gen<CuC,double,N>(d_res, 1.0, d_mtmp, d_res);`.
- **`adj` (`:405`)**: identical diagonal add ($M$ real $\Rightarrow m_L^\dagger=m_L$).
- **squared op (`:462-503`)** and **`_ms` variant (`:494-503`)**: implement
  $(1+M)(D_\text{ov}v)+D_\text{ov}^\dagger[(1+M)v]+M^2v$ using `M_mass`/`mass_coeff` for the
  $M\cdot$ applies (the old code already computes `mult`/`adj`; reuse them).
- **mass inner products (`:633-692`)**: replace the scalar-`mass` combinations consistently
  (only needed if the L=1 test also exercises the action/force; for a pure mult/adj/Dsq
  check this can come later).

### 3. Test program `test_diag_mass_l1_claude.cu`
Mirror the L=1 driver setup (`hmc_fermilab_wmass_L2_claude.cu` with `N_REFINE=1`):
build `Base base(1); WilsonDirac DW(base, 0.0, 1.0, M5, at, nu0);`. Then:
```
#include "includes/overlap_wmass_obsolete_claude.h"   // obsolete::OverlapWMass (scalar)
#include "includes/overlap_wmass_claude.h"            // OverlapWMass (diagonal)
...
const double m_phys = 0.1;                            // any test value
const double c = m_phys * base.mean_dual_area / base.mean_ell;   // uniform L=1 value
obsolete::OverlapWMass<WilsonDirac> Dold(DW, Complex(c,0.0), npole);  Dold.update(U);
OverlapWMass<WilsonDirac>           Dnew(DW, Complex(m_phys,0.0), npole); Dnew.update(U);
// random vector v (fixed seed); apply and compare:
//   mult:  ||Dnew.mult(v)  - Dold.mult(v)||  / ||Dold.mult(v)||
//   adj:   same with adj
//   Dsq:   same with the squared-operator apply
```
Use a free gauge (`U` identity / as the eig test) — the mass reduction is gauge-independent,
identity is simplest.

### 4. Build & run (local)
Add a `Makefile.fnal` (or local makefile) target for `test_diag_mass_l1_claude.o`
mirroring the hmc targets (`N_REFINE=1`). Run on a GPU. **Expected:** all three relative
differences $\sim 10^{-13}$–$10^{-14}$ (machine precision); `mass_diag_h`/`M_mass` uniform
at $L=1$.

## Step 3 (expanded scope, 2026-06-19)

**(a) Operator equivalence — BOTH paths.** For each m in {0, 0.1, 0.1 i}, compare
obsolete(scalar c) vs production(physical m) on a Gaussian source for the **non-ms**
`mult`/`adj`/`DHD` AND the **`_ms`** `mult_ms`/`adj_ms`/`DHD_ms`. (Both must reduce to the
scalar reference at L=1 to ~1e-13.)

**(b) HMC force vs numerical derivative** (mimics the massive-operator check in
`hmc_fermilab_claude.cu:357-487`; now in `test_diag_mass_l1_claude.cu`). The mass enters
the force as the $(1+M^*)$ weight on the force inner products (`grad_*`), so validate directly:
- Fix a **Gaussian pseudofermion** $\phi$ (set `d_phi` once; do NOT re-heatbath while
  varying U, since the heatbath uses $D_m$). Action $S(U)=\phi^\dagger(D_m^\dagger D_m)^{-1}\phi$
  = `update_eta()` then `action()`.
- Analytic per-link force: `D.grad_deviceAsyncLaunch(link, U, d_eta)` (with `d_eta` from the
  unperturbed U). Pick **a few links** (e.g. 3-4 spatial + 1 temporal).
- Numerical: bump that link's gauge angle $\theta_\ell\to\theta_\ell\pm\epsilon$ (eps~1e-4),
  recompute $S$, central difference $\partial S/\partial\theta_\ell\approx[S(+)-S(-)]/2\epsilon$.
  Compare to the analytic grad (mind the force sign convention = $-\partial S/\partial\theta$).
- Do this for the **production** operator (validates the diagonal grad) AND check it equals
  the **obsolete** force at L=1. Tolerance ~1e-5..1e-6 (FD truncation-limited).
- PREREQUISITE: the `grad_*` routines must be converted first (fold $\bar M$ into the dots;
  the active variant is `_l4` via `#define GRAD_L4`; the `_l2`/`_l4` block-dot fold needs a
  block $\bar M$ weight). For the FD check, simplest is the **default** `grad` (build WITHOUT
  `GRAD_L1/2/4`) once it is converted.

## Pass criterion
mult, adj, and $D_m^\dagger D_m$ all agree obsolete-vs-production to $\sim10^{-14}$ at $L=1$.
On pass → proceed to the driver/CLI chunk (pass physical $m$; reverse-engineer the four
production masses $m=m_1\,\bar a_s^{(1)}/A_y^{(1)}$) and an $L=2$ sanity run.

## Remaining work & difficulties (which routine, why)

**DONE:** `mult`/`adj`/`DHD` (+ `_ms`), and the **default** `grad_deviceAsyncLaunch`
(per-pole streamed) -- folded via a stream-aware `M_mass` apply (`MatPoly Mp` bound to
`stream[istream]`) into `d_Ms[m]` + an extra `dotAsync`, with `conj(mass_coeff)` folded
host-side. No device-sync inside the omp pole loop.

**TODO -- block force variants `grad_*_l1/_l2/_l4`** (`overlap_wmass_claude.h`, ~lines
770-916). The ACTIVE production force is **`_l4`** (`#define GRAD_L4`), so this gates HMC.
- `grad_l1`: **easy.** Same per-pole streamed structure as the default grad (just reads
  precomputed `d_XZpre/d_XYpre`); apply the identical stream-aware `M_mass` dot fold.
- `grad_l2` / `grad_l4`: **harder -- this is the real work.** They have NO per-pole loop:
  the poles are a contiguous `N*npole` block (`d_XZg/d_XYg/d_CY/d_CZ`) and the dots are one
  batched `block_dot`. To fold `(1+M^*)` you must:
  1. write a **block diagonal-weight kernel** that multiplies the block element
     `[m*N+i] *= (A_y/Abar)_i` (i.e. `M_mass` broadcast across all `npole` columns) -- apply
     it to the X-applied block (`d_XZg`/`d_XYg`), into a scratch block;
  2. run a **second `block_dot`** of `d_CY`/`d_CZ` against that weighted block -> `d_dotA_M`/`d_dotB_M`;
  3. in the host reduce, use `real(am + conj(mass_coeff)*am_M)` (and same for `bm`).
  Plus the post-loop m=0 term needs the same fold (as in the default grad).
  Difficulty = the new block-diagonal kernel matching `block_dot`'s `[m*N+i]` layout +
  two extra `npole`-length dot buffers; `_l4` also routes the link via `link_matvec_block`
  but the M-weight still goes on the `XZg/XYg` (X-applied) blocks, not `CY/CZ`.

**TODO -- driver/CLI** `hmc_fermilab_wmass_L{2,4}_claude.cu`: the `mass_re/mass_im` CLI arg
is now the **physical** m (R=1 units), not the old bare value. Recompute the production
masses `m = m_1 * abar_s^(1)/A_y^(1)` and update the wrapper `PAIR_LIGHT/HEAVY`. Gotcha:
`dir3` (checkpoint path) encodes `mass.real()/imag()`, so the physical-m runs land in NEW
checkpoint dirs (new ensembles) -- the old bare-mass dirs are not reused.

## Notes / gotchas
- Namespace-isolate the obsolete header or the two `OverlapWMass`/`Zolotarev` definitions
  collide at link/compile.
- `mass_coeff` uses `mean_dual_area` (not `mean_ell`) in the numerator ratio — it converts
  `volume_matrix`'s $A_y/\bar A$ to the target $A_y/\bar a_s$.
- Keep `mass` real in the test; the complex-conjugate paths (`adj` uses `conj(mass)`) are
  trivial for real $M$ but keep them correct for the general case.
- This is the gauge-independent operator check only; the HMC force is unchanged (decision 6).
