# Plan: Connected and Disconnected Current Correlator Measurements

## Physics / Goal Summary

This plan specifies two new measurement programs that compute current-current two-point
functions using the exactly conserved overlap vector current $J_{V,+}^{nn'}(t)$, defined
on the link $(n,n')$ at time $t$.  The correlators are formed by Wick contraction of two
current insertions; for the $N_f = 2$ case they split into a connected diagram ($C_{V,c}$)
and a disconnected diagram ($C_{V,d}$).

The derivation and notation follow Section 3.3 of
`/mnt/barracuda22/qed3/qed3int_v2.pdf` (Katz-Matsumoto, dated 2026-06-01), Eqs. (3.39)--(3.55).

---

## Physics Background

### Notation (two-component, $N_f = 2$, $\xi$/$\eta$ basis)

The Wick contraction rules (Eqs. (3.39)--(3.40)) are

$$\overline{\xi_{f,x_1}\eta_{g,x_2}^\dagger} = \delta_{fg}(D_\text{ov}^{-1})_{x_1 x_2}, \qquad
  \overline{\eta_{f,x_1}\xi_{g,x_2}^\dagger} = \delta_{fg}(D_\text{ov}^{-\dagger})_{x_1 x_2}.$$

For a spatial link $(n,n')$ on timeslice $t$, the kernel is $K^{nn'} \equiv K^{wz}$
(two-component, Eq. (3.32) and its adjoint), implemented in
`includes/conserved_current_claude.h::ConservedCurrent::apply_k`.

The (+) branch of the vector current (two-component, Eq. (3.36)) is

$$J_{V,+}^{nn'}(t) = \eta_f^\dagger K^{nn'}(t) \xi_f.$$

The current correlators at time separation $\Delta t$ summed over source timeslice $t_0$
and over all links $(n,n')$ within the same face (see Eq. (3.55)) are defined below.

### Connected vector correlator $C_{V,c}^{nn'}(0 \to t; U)$ (Eq. (3.41))

$$
\frac{2}{N_f}\langle J_{V,+}^{nn'}(0)\,J_{V,+}^{nn'}(t)\rangle_F
  = -\mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\,D_\text{ov}^{-1}K^{nn'}(t)\right]
  + \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\right]\mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(t)\right]
$$

$$\equiv -C_{V,c}^{nn'}(0\to t;U) + C_{V,d}^{nn'}(0\to t;U). \tag{3.41, 3.42}$$

The connected piece is

$$C_{V,c}^{nn'}(0\to t;U) = \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\,D_\text{ov}^{-1}K^{nn'}(t)\right]. \tag{3.41}$$

The disconnected piece is

$$C_{V,d}^{nn'}(0\to t;U)
  = \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\right]
    \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(t)\right]. \tag{3.42}$$

Note: $\mathrm{tr}[D_\text{ov}^{-1}K^{nn'}] = -J_{V,+}^{nn'}/(N_f)$ from the formula
$J_{V,+} = -N_f\,\mathrm{tr}(K D_\text{ov}^{-1})$ already verified in
`check_conserved_current_claude.cu` Chunks 5b and 6.

### Symmetry under time reversal (Eq. (3.43)--(3.44))

$$
\frac{2}{N_f}\langle J_{V,-}^{nn'}(0)\,J_{V,-}^{nn'}(t)\rangle_F
  = -\bigl(C_{V,c}^{nn'}(0\to t;U)\bigr)^* + \bigl(C_{V,d}^{nn'}(0\to t;U)\bigr)^*. \tag{3.43, 3.44}
$$

This gives a consistency check: $C_{V,c}(t)^* = C_{V,c}(N_t - t)$ etc.

### Axial correlator $C_A^{nn'}$ (Eqs. (3.45)--(3.48))

$$\frac{2}{N_f}\langle J_{A,+}^{nn'}(0)\,J_{A,-}^{nn'}(t)\rangle_F
  = \mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)(1-D_\text{ov}^\dagger)\,
      D_\text{ov}^{-\dagger}K^{nn'\dagger}(t)(1-D_\text{ov})\right]
  \equiv C_A^{nn'}(0\to t;U). \tag{3.45, 3.46}$$

The axial correlator involves $(1 - D_\text{ov})$ insertions and is purely connected
(no disconnected term from the Wick contraction).

### Flavor current correlator (Eqs. (3.49)--(3.52))

For the flavor $SU(N_f/2)^2$ current (Eq. (3.38), (+) branch):

$$\langle J_{r,+}^{nn'}(0)\,J_{r',+}^{nn'}(t)\rangle_F
  = -\delta_{rr'}\,\mathrm{tr}\!\left[D_\text{ov}^{-1}K^{nn'}(0)\,D_\text{ov}^{-1}K^{nn'}(t)\right]
  = -\delta_{rr'}\,C_{V,c}^{nn'}(0\to t;U). \tag{3.49, 3.50}$$

The flavor correlator is identical to the connected vector piece; no separate computation needed.

### Stochastic estimators

Both traces $\mathrm{tr}[D^{-1}K^{nn'}\cdot D^{-1}K^{nn'}]$ (connected) and
$\mathrm{tr}[D^{-1}K^{nn'}]$ (disconnected loop) are estimated with $Z_2 \times Z_2$ noise.

For the **disconnected loop** at link $(n,n')$, timeslice $t$:

$$L^{nn'}(t) \equiv -N_f\,\mathrm{tr}(K^{nn'}(t)\,D_\text{ov}^{-1})
  \approx -N_f \cdot \frac{1}{N_\text{hits}} \sum_h \eta_h^\dagger K^{nn'}(t)\,\phi_h,$$

where $D_\text{ov}\phi_h = \eta_h$ (solved via normal equations $D^\dagger D\phi = D^\dagger\eta$).
This is exactly the estimator already working in `check_conserved_current_claude.cu` Chunk 5b.

For the **connected piece**, the one-end trick gives:

$$C_{V,c}^{nn'}(0\to t) \approx \frac{1}{N_\text{hits}}\sum_h
  \bigl(\eta_h^\dagger K^{nn'}(0)\,\phi_h\bigr)
  \cdot \overline{\bigl(\eta_h^\dagger K^{nn'}(t)\,\phi_h\bigr)},$$

where the same $\phi_h$ from the same solve $D\phi_h = \eta_h$ is reused.
This is an unbiased estimator of $\mathrm{tr}[D^{-1}K^{nn'}(0)\cdot D^{-1}K^{nn'}(t)]$
via $E[\eta^\dagger A\phi \cdot \overline{\eta^\dagger B\phi}] = \mathrm{tr}[A D^{-1} (B D^{-1})^\dagger]$
for $Z_2$ noise.

**Note on the $(0 \to t)$ structure**: source is at $t_0 = 0$ (wall), sink at $t$. In
practice we average over source timeslice $t_0$ (see Eq. (3.55)) and wrap modulo $N_t$.

### Lattice-to-continuum normalization (Eqs. (3.54)--(3.55))

The renormalization-free lattice-to-continuum translation uses the simplicial geometry:

$$\frac{1}{A_\triangle} \sum_{(n,n') \in \triangle}
    \frac{\ell_{nn'}^3 \ell_{nn'}^*}{p_0^2 A_{nn'}^2}
    J_V^{nn'}(0) J_V^{nn'}(t)
    \sim (\delta^{ab} - e_3^a e_3^b) j_V^a(x_\triangle, 0)\, j_V^b(x_\triangle, t). \tag{3.55}$$

This geometric prefactor is applied at the analysis stage in Julia/Python; the measurement
program outputs un-normalized per-link correlator values.

---

## References

- A. Katz, N. Matsumoto, "Notes on QED3 (interacting theory)", dated 2026-06-01,
  `/mnt/barracuda22/qed3/qed3int_v2.pdf` --- Sec. 3.3 (current correlators on the lattice,
  Eqs. 3.39--3.55); Sec. 3.2.2 (two-component formalism for $K^{wz}$, Eq. 3.32).
  This is the primary physics reference for this plan.

---

## Scope and Constraints

**New files to create:**

| File | Role |
|---|---|
| `src/both_3d/jj_conn_claude.cu` | Connected vector correlator $C_{V,c}$ measurement driver |
| `src/both_3d/jj_disc_claude.cu` | Disconnected vector correlator loop $L^{nn'}(t)$ measurement driver |

**Existing files used (not modified):**

| File | Role |
|---|---|
| `includes/conserved_current_claude.h` | `ConservedCurrent::apply_k` --- the kernel $K^{wz}$ already implemented and verified |
| `includes/valence_claude.h` | `FermionVector`, `fill_z2_source`, `dag`, `time_spin_dilution`, `accumulate_loop_gamma` |
| `includes/overlap.h` | `Overlap<WilsonDirac>`, `adj_deviceAsyncLaunch`, `DHD_deviceAsyncLaunch` |
| `disc_claude.cu` | Template for disconnected loop structure, HDF5 output, CG solve pattern |
| `meson_pq_wall_v2_claude.cu` | Template for connected correlator outer loop, HDF5 output |
| `check_conserved_current_claude.cu` | Template for CG solve chain with `ConservedCurrent` |

**Compile-time parameters (both programs inherit from `disc_claude.cu` defaults):**

```cpp
namespace Comp {
  constexpr int N_REFINE = 1;
  constexpr int Nt       = 128;
  constexpr int NPARALLEL_DUPDATE = 1;  // set to 1; increase for production
  const double TOL_INNER = 1.0e-15;
  const double TOL_OUTER = 1.0e-14;
}
```

**Runtime parameters (getopt_long, mirroring `disc_claude.cu`):**

| Flag | Default | Meaning |
|---|---|---|
| `--gsq` | 8.0 | Wilson coupling $g^2$ |
| `--Nf` | 2 | Number of sea quark flavors |
| `--nu0` | 1.0 | Sea quark asymmetry $\nu_0$ |
| `--nu1` | 1.0 | Valence quark asymmetry $\nu_1$ |
| `--nhits` | 1 | Number of $Z_2\times Z_2$ noise hits per config |
| `--t_block` | 8 | Timeslices per dilution block (disconnected only) |

**Out of scope:**

- Axial correlator $C_A$ (Eqs. 3.45--3.48): deferred; requires $(1-D_\text{ov})$ insertions.
- Flavor correlator $C_{r}$: identical to $C_{V,c}$, no separate program needed.
- Temporal links in the correlator: implementation uses spatial links only (same scope
  as the existing `apply_k_wz` wrapper).
- Julia analysis scripts for the geometric normalization (Eq. 3.55): out of scope here.
- $L=2$ (N_REFINE=2) variant: out of scope; same code compiles with different `N_REFINE`.

---

## Data Structures

### Reuse: `FermionVector` (`includes/valence_claude.h`)

Host-pinned, `field[N]` array of `Complex`.  Methods used:

- `fill_z2_source(rng)` --- all-to-all $Z_2\times Z_2$ noise (connected program)
- `time_spin_dilution(rng, t_s, t_block, spin)` --- diluted source (disconnected program,
  mirroring `disc_claude.cu:293`)
- `dag(b)` --- host inner product $\sum_i \bar{a}_i b_i$ (from `valence_claude.h:~186`)

### Reuse: `ConservedCurrent<Fermion>` (`includes/conserved_current_claude.h`)

Already constructed and verified.  Key method:

```cpp
// TEMPLATE: check_conserved_current_claude.cu:879
kop.apply_k(d_k_xi, d_phi_dev, U, std::pair<int,BaseLink>{s_focus, lk_j});
CUDA_CHECK(cudaDeviceSynchronize());
CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi.field), d_k_xi, N*CD, D2H));
const Complex c = eta.dag(fv_Kphi);
```

### New per-program: per-link accumulator arrays

For `jj_disc_claude.cu`: `std::vector<std::vector<Complex>> L(n_links, std::vector<Complex>(Nt, 0.0))`
accumulates the per-hit stochastic estimate of $-N_f \mathrm{tr}(K^{nn'}(t) D^{-1})$.

For `jj_conn_claude.cu`: `std::vector<std::vector<Complex>> Cc(n_links, std::vector<Complex>(Nt, 0.0))`
accumulates the per-hit estimate of $\mathrm{tr}[D^{-1}K^{nn'}(0)\cdot D^{-1}K^{nn'}(t)]$.

---

## Ordered Implementation Chunks

### Chunk 1 --- `jj_disc_claude.cu`: boilerplate and construction

**What:** Copy the boilerplate of `disc_claude.cu` (includes, Comp namespace, type aliases,
CLI parsing, GPU setup, lattice/DW/Overlap construction, `ConservedCurrent` construction,
`d_MemorySets` allocation).  The only addition is `#include "includes/conserved_current_claude.h"`.

**Key function to write:** `ParseArgs` --- identical to `disc_claude.cu:105-139`.

Template:

```cpp
// TEMPLATE: disc_claude.cu:1-88 (includes, Comp namespace, type aliases)
// TEMPLATE: disc_claude.cu:156-229 (main: GPU setup, Base, DW, Fermion, Gauge construction)
// TEMPLATE: check_conserved_current_claude.cu:226-228 (ConservedCurrent construction)
ConservedCurrent<Fermion> kop(D);
```

**Files:** `jj_disc_claude.cu`

---

### Chunk 2 --- `jj_disc_claude.cu`: enumerate links and define output layout

**What:** Build `link_list` --- a `std::vector<std::pair<int,BaseLink>>` enumerating all
spatial links across all timeslices, i.e.

```cpp
for(int s = 0; s < Nt; s++)
  for(const auto& lk : base.links)
    link_list.push_back({s, lk});
```

This gives `n_links = Nt * base.n_links` entries.

Define output directory `dir4` and HDF5 file path following `disc_claude.cu:206-214`.
Suggested path:

```
jj_disc_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nu1{nu1}nt{Nt}L{N_REFINE}tb{t_block}/jj_disc_corr.{k}.h5
```

HDF5 dataset layout: one dataset per link, keyed by `"{h}/{s}/{ix0}/{ix1}/real"` and
`".../imag"`, each of length `Nt`.

Resume sentinel (mirrors `disc_claude.cu:272-278`): check existence of the last expected
dataset before processing each $k$.

**Files:** `jj_disc_claude.cu`

Template:

```cpp
// TEMPLATE: disc_claude.cu:204-278 (dir3/dir4 construction, HDF5 resume logic)
```

---

### Chunk 3 --- `jj_disc_claude.cu`: per-config stochastic loop

**What:** For each gauge config $k$, for each hit $h$:

1. Draw diluted source `eta.time_spin_dilution(rng, t_s, t_block, spin)`.
2. Solve $D\phi = \eta$ via normal equations (identical to `disc_claude.cu:294-295`):

```cpp
// TEMPLATE: disc_claude.cu:294-295
op_DH.from_cpu<N>(DH_eta.field, eta.field);
op_DHD.solve<N>(phi.field, DH_eta.field);
CUDA_CHECK(cudaMemcpy(d_phi_dev, reinterpret_cast<const CuC*>(phi.field), N*CD, H2D));
```

3. For each spatial link $(s, lk)$ in `link_list`:

```cpp
kop.apply_k(d_k_xi, d_phi_dev, U, std::pair<int,BaseLink>{s, lk});
CUDA_CHECK(cudaDeviceSynchronize());
CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi.field), d_k_xi, N*CD, D2H));
const Complex c = eta.dag(fv_Kphi);
L_acc[link_idx][s] += c;  // accumulate over dilution/hits
```

4. After all $(t_s, \text{spin})$ dilution blocks and all hits, apply the $-N_f$ factor and
   compute the disconnected correlator $G_\text{disc}[\Delta t] = (1/N_t)\sum_{t_0} L(t_0) L^*(t_0+\Delta t)$
   using `compute_disc_corr` (same as `disc_claude.cu:142-154`).

5. Write per-link real and imaginary parts to HDF5.

**Efficiency note:** `apply_k` for each link involves `D.size-1` CG solves (one per
Zolotarev pole).  For `N_REFINE=1` (`base.n_links=30`) and `Nt=128`, there are $30\times 128 = 3840$
link calls per solve.  This will be slow; see Open Questions on batching.

**Files:** `jj_disc_claude.cu`

Template:

```cpp
// TEMPLATE: disc_claude.cu:288-308 (hit loop, dilution, CG solve, accumulate)
// TEMPLATE: check_conserved_current_claude.cu:870-885 (apply_k, dag inner product)
```

---

### Chunk 4 --- `jj_disc_claude.cu`: HDF5 output and cleanup

**What:** After all hits for config $k$, write output to HDF5.  Dataset keys:
`"{h}/{s}/{lk[0]}/{lk[1]}/real"` and `".../imag"`, each a `std::vector<double>` of length `Nt`.
Pattern mirrors `disc_claude.cu:303-309`:

```cpp
// TEMPLATE: disc_claude.cu:303-309
h5file.createDataSet(ds_prefix + "real", G_real[link_idx]);
h5file.createDataSet(ds_prefix + "imag", G_imag[link_idx]);
```

Cleanup: `cudaFree(d_phi_dev)`, `d_MemorySets[i].deallocate()`.

**Files:** `jj_disc_claude.cu`

---

### Chunk 5 --- `jj_conn_claude.cu`: boilerplate and construction

**What:** Same structure as Chunk 1 for the connected program.  Identical includes, Comp
namespace, CLI, GPU setup, and operator construction.  No `--t_block` parameter (connected
uses all-to-all noise, no time dilution).

**Files:** `jj_conn_claude.cu`

Template:

```cpp
// TEMPLATE: disc_claude.cu:1-229 (full boilerplate)
// TEMPLATE: check_conserved_current_claude.cu:226-228 (ConservedCurrent construction)
```

---

### Chunk 6 --- `jj_conn_claude.cu`: per-config stochastic loop (one-end trick)

**What:** For each gauge config $k$, for each hit $h$:

1. Draw noise `eta.fill_z2_source(rng)` (all-to-all, no dilution).
2. Solve $D\phi = \eta$ once per hit via normal equations.
3. For each link $(s, lk)$:

```cpp
kop.apply_k(d_k_xi, d_phi_dev, U, std::pair<int,BaseLink>{s, lk});
CUDA_CHECK(cudaDeviceSynchronize());
CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(fv_Kphi.field), d_k_xi, N*CD, D2H));
const Complex c0 = eta.dag(fv_Kphi);   // estimator of tr(K^{nn'} D^{-1}) at time s
```

4. The connected correlator between timeslices $s$ and $(s + \Delta t) \bmod N_t$:

$$C_{V,c}^{nn'}(\Delta t) \approx c_0(s) \cdot \overline{c_0((s+\Delta t)\bmod N_t)}.$$

This requires storing `c0[s]` for all $s$ for each link, then computing the outer product
in $\Delta t$.  For each link, accumulate `Cc_acc[link_idx][dt] += c0[s] * conj(c0[(s+dt)%Nt])`.

5. After all hits, divide by `nhits * Nt`, apply factor $-N_f$, and write.

**Efficiency note:** storing `c0[Nt]` per link per hit requires `n_links * Nt` complex
scalars in host memory; for `N_REFINE=1`, that is `30*128 * 128 = 491520` complex doubles
$\approx$ 7.5 MB --- negligible.

**Files:** `jj_conn_claude.cu`

Template:

```cpp
// TEMPLATE: check_conserved_current_claude.cu:870-885 (apply_k, dag, per-link accumulation)
// TEMPLATE: disc_claude.cu:142-154 (compute_disc_corr --- same structure for C_c)
```

---

### Chunk 7 --- `jj_conn_claude.cu`: HDF5 output

**What:** Write per-link correlator to HDF5.  Suggested dataset key:
`"{h}/{s}/{lk[0]}/{lk[1]}/real"` and `".../imag"`, each length `Nt`.
Output directory:

```
jj_conn_Nf{Nf}_gsq{gsq}at{at}nu0{nu0}nu1{nu1}nt{Nt}L{N_REFINE}/jj_conn_corr.{k}.h5
```

Resume logic same as Chunks 2/4.

**Files:** `jj_conn_claude.cu`

---

### Chunk 8 --- Validation check

**What:** Add a small self-consistency check: for `N_REFINE=1`, `Nt=4`, compare the
stochastic $C_{V,d}^{nn'}(t)$ from `jj_disc_claude.cu` against the exact basis-vector
computation (same guard `if constexpr (N_REFINE==1 && Nt<=4)` as in
`check_conserved_current_claude.cu:643`).  For the connected piece, verify that
$\sum_{nn'} C_{V,c}^{nn'}(0) \geq 0$ (positivity at zero separation).

This can be implemented as a compile-time-guarded block at the end of the `jj_conn_claude.cu`
main, printing `# PASS / FAIL` lines, before the production `Nt=128` run.

**Files:** `jj_conn_claude.cu` (guarded block only; no new file)

---

## Interface to Existing Code

### `ConservedCurrent::apply_k` (already implemented, not modified)

The generic template `apply_k<Gauge, LinkEl>` in `includes/conserved_current_claude.h`
accepts `std::pair<int,BaseLink>` for spatial links.  The calling pattern established in
`check_conserved_current_claude.cu` is the canonical interface; both new programs follow it.

### `MatPoly::solve` / `from_cpu` (not modified)

CG solve chain:

```cpp
// TEMPLATE: disc_claude.cu:231-237
auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch,  &D, ...);
auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch,  &D, ...);
LinOpWrapper M_DH(f_DH), M_DHD(f_DHD);
MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
MatPoly op_DHD; op_DHD.push_back(cplx(1.0), {&M_DHD});
// ...
op_DH.from_cpu<N>(DH_eta.field, eta.field);
op_DHD.solve<N>(phi.field, DH_eta.field);
```

### `GaugeExt::read` (not modified)

```cpp
// TEMPLATE: disc_claude.cu:282
U.read(dir3 + "ckpoint_lat." + std::to_string(k));
D.update(U);
```

### `FermionVector::dag` (not modified)

```cpp
// TEMPLATE: check_conserved_current_claude.cu:259
const Complex c = eta.dag(fv_Kphi);
```

---

## File Naming

All new files follow the `_claude` suffix convention:

- `src/both_3d/jj_disc_claude.cu`
- `src/both_3d/jj_conn_claude.cu`

---

## Open Questions

1. **Link averaging vs per-link output.** The correlators depend on the link index $(n,n')$.
   For physical observables, one projects onto irreps of the icosahedral symmetry group
   (orbit-averaged).  Should the programs output:
   (a) per-link data (maximum flexibility, larger HDF5 files), or
   (b) orbit-averaged correlators (smaller output, but requires knowing the orbit decomposition)?
   Currently planning (a).

2. **Number of links per output.** For `N_REFINE=1`, `base.n_links=30`, `Nt=128`, there are
   3840 spatial links total.  Storing `C_c(Nt)` per link per config produces
   $3840 \times 128 \times 2 \times 8 = 7.5$ MB per config.  Is this acceptable, or should
   we output only a subset of links (e.g., one per orbit)?

3. **Dilution for the connected piece.** The current plan uses all-to-all noise for
   `jj_conn_claude.cu` (no time dilution).  For better signal at large $\Delta t$,
   time-spin dilution (as in `disc_claude.cu`) could be applied.  Confirm whether
   all-to-all is sufficient or dilution is needed.

4. **Temporal current contributions.** Eqs. (3.41)--(3.52) are stated for spatial links
   $(n,n')$; `apply_k_wz` (wrapped by `apply_k`) handles temporal links via the
   `std::pair<int,Idx>` overload.  Should the temporal-link correlators also be measured,
   or spatial-only for now?

5. **Source timeslice averaging.** The plan averages over source timeslice $t_0$ by
   wrapping indices mod $N_t$ (same as `compute_disc_corr`).  Confirm this is the correct
   interpretation of Eq. (3.55), which sums face-projected correlators.

6. **Axial correlator.** Eqs. (3.45)--(3.48) give $C_A^{nn'}$ which requires
   $(1 - D_\text{ov})$ applied to the right of $D_\text{ov}^{-1}K^{nn'}$.  This is an
   additional operator application per pole.  Implement now alongside $C_{V,c}$, or defer?

7. **Production `nhits` and dilution scheme.** For `disc_claude.cu`, `t_block=8` and
   `nhits=1` is the production setting.  What are the target `nhits` and noise scheme for
   the new programs?  This affects expected runtime estimates.

---

## Efficiency Decisions Recorded

- **One-end trick for the connected piece (Chunk 6):** reuse the single solve $D\phi=\eta$
  for both the source-end ($s=0$) and the sink-end ($s=\Delta t$) kernels.  This halves
  the number of solves vs a two-end approach.  Cost: $N_\text{hits}$ solves per config.

- **Per-link sequential `apply_k` calls:** each call to `apply_k` for one link launches
  `D.size-1` CG solves sequentially.  For 3840 links and 21 poles, the total CG calls
  per hit per config is $3840 \times 20 = 76800$.  This is the dominant runtime cost.
  Batching over links (reusing pole-$\ell$ solutions across multiple links on the same
  timeslice) would reduce this, but is left as a future optimization.  Decision: use
  the simple sequential approach first; profile before optimizing.

- **No time dilution for connected piece:** all-to-all $Z_2\times Z_2$ noise is used in
  `jj_conn_claude.cu`; the timeslice structure of $c_0(s)$ is extracted from the single
  propagator $\phi$ per hit, not from separate per-timeslice solves.

- **HDF5 per-trajectory, Truncate mode:** one `.h5` file per trajectory checkpoint $k$,
  opened with `HighFive::File::Truncate`, mirroring `disc_claude.cu:285-286`.
