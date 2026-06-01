# Implementation Plan: Check Program for the Exactly Conserved Vector Current $K^{wz}$

_Created: 2026-05-25_

---

## 1. Physics Motivation

The kernel `apply_k_wz` in `conserved_current_claude.h` implements $K^{wz}$, the
two-component overlap-fermion exactly conserved vector current for a single spatial
link $(s \to \text{lk})$. Its correctness is verified via:

1. **Theta checks (passed).** $[X,\Theta]\xi = \sum_l \delta\theta_l W^l\xi$ and
   $[D_\text{ov},\Theta]\xi = \sum_l \delta\theta_l K^l\xi$ hold to machine precision for
   random $\Theta$, confirming $K^{wz}$ is the exact link derivative of $D_\text{ov}$.

2. **Stochastic Ward identity (Chunk 5b).** $\sum_{z:\text{nn of }w}\langle J_V^{wz}\rangle = 0$
   verified stochastically via Z2$\times$Z2 noise.

Note: $K^{wz}$ is NOT Hermitian (see Resolved Decisions).

---

## 2. Files

### Created (new)
- `src/both_3d/check_conserved_current_claude.cu` — standalone check program

### Read-only references (not modified)
- `src/both_3d/includes/conserved_current_claude.h` — kernel under test
- `src/both_3d/includes/overlap.h` — `Overlap<WilsonDirac>` type
- `src/both_3d/includes/dirac_simp.h` — `DiracS2Simp`, `d_coo_format`
- `src/both_3d/includes/gauge_ext.h` — `GaugeExt`, `read()`, `base.links`
- `src/both_3d/includes/sparse_matrix.h` — `COO<N>`
- `src/both_3d/includes/gpu_header.h` — `Taxpy_gen`, `NBlocks`, `CD`, `D2H`, `H2D`
- `src/both_3d/disc_claude.cu` — config-loading and Overlap setup reference

---

## 3. Comp Parameters

Mirror `disc_claude.cu:36--57`. Set `NSTREAMS=1` for single-stream diagnostics:

```cpp
// In check_conserved_current_claude.cu
namespace Comp {
  constexpr int N_REFINE    = 1;
  constexpr int Nt          = 128;
  constexpr int NS          = 2;
  // N_SITES, N derived from N_REFINE (10*N_REFINE^2 + 2 = 12 sites)
  const double TOL_INNER    = 1.0e-9;
  const double TOL_OUTER    = 1.0e-8;
}
```

---

## 4. Ordered Implementation Chunks

---

### Chunk 1 — Boilerplate: includes, Comp namespace, type aliases, main() skeleton

Follow `disc_claude.cu:1--183`. Differences:
- No HDF5 / HighFive include.
- Add `#include "includes/conserved_current_claude.h"` after the overlap include.
- Command-line args: `gsq`, `Nf`, `nu0` (identical positional order as `disc_claude.cu`).
- Mirror `d_MemorySets` allocation if required by `sparse_dirac.h` (see Open Question 3).

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 2 — Lattice, WilsonDirac, Overlap, and ConservedCurrent construction

```cpp
// (inside main(), after GPU setup)
using Base         = S2Simp;
using WilsonDirac  = DiracExt<Base, DiracS2Simp>;
using Gauge        = GaugeExt<Base, Comp::Nt, Comp::is_compact>;
using Fermion      = Overlap<WilsonDirac>;

Base         base(Comp::N_REFINE);
const double at   = 0.2;
WilsonDirac  DW(base, 0.0, 1.0, /*M5=*/-1.0, at, nu0);
Fermion      D(DW, 21);
Gauge        U(base);       // default-constructed: free-field (zero phases)

D.update(U);

ConservedCurrent<Fermion> kop(D);
```

Allocate device vectors `d_xi`, `d_k_xi` of size `Comp::N * CD`.

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 3 — Helper functions

Three host-side helpers:

```cpp
// inner product <eta, v> = sum_i conj(eta[i]) * v[i]
Complex inner_product(const std::vector<Complex>& eta,
                      const std::vector<Complex>& v);

// transfer device CuC buffer to host Complex vector of length N
void device_to_host_vec(std::vector<Complex>& h_res,
                        const CuC* d_res, int N);

// accumulate lattice divergence for one link (s_lat, lk)
// h_k_xi: host copy of K^{wz}\xi for this link, length Comp::N
// sign convention: already embedded in apply_k_wz output via W^{wz} antisymmetry
// (see dirac_simp.h:294--338; plain sum gives the lattice divergence)
void accumulate_divergence(std::vector<Complex>& h_div,
                           const std::vector<Complex>& h_k_xi,
                           int s_lat, const BaseLink& lk);
```

The sign in `accumulate_divergence` must be verified against the COO entries produced
by `build_W_wz` before finalizing (see Open Question 1).

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 4 — Cross-check via HMC force code

**Mathematical basis.** Working through `grad_deviceAsyncLaunch` (overlap.h:423--474)
term by term, using `d_coo_format` = $\partial D_W/\partial\theta_{wz}$ = $I W^{wz}$ and
$\text{Re}[i\langle a, W^{wz}b\rangle] = -\text{Im}[\langle a, W^{wz}b\rangle]$:

$$F^{wz}(\eta) = \frac{2}{\lambda_\text{max}} \text{Im}[\langle\eta, K^{wz}\eta\rangle]$$

This identity holds term-by-term for both Term A and each pole in the sum.
The test compares `Im(<xi, K^{wz} xi>)` from `apply_k_wz` against
`lambda_max/2 * grad_deviceAsyncLaunch` for every link.

**Why this is useful.** `grad_deviceAsyncLaunch` is a completely different code path
(used in HMC, independently written), so this provides genuine cross-validation.
It tests the imaginary part of the diagonal matrix element per link; only a scalar
check (not a vector check), but independent of Chunks 5 and 6.

**Procedure.**

1. Call `D.precalc_grad_deviceAsyncLaunch(U, d_xi)` once (free-field $U$).
   This sets `D.d_Zs[m]` = $R_m\xi$, `D.d_Ys[m]` = $R_m X^\dagger\xi$ for $m = 1\ldots\text{size}-1$.
2. For `s_lat` in $[0, N_t)$, for `lk` in `base.links`:
   - `force = D.grad_deviceAsyncLaunch(std::pair<int,BaseLink>{s_lat,lk}, U, d_xi)`;
   - `kop.apply_k_wz(d_k_xi, d_xi, U, s_lat, lk)`; sync; copy to host;
   - `inner = inner_product(h_xi, h_k_xi)`;
   - `err = |Im(inner) - 0.5 * D.lambda_max * force|`.
3. Report `max_err` over all links; assert `< 1.0e-4`.

**Note on scratch.** `apply_k_wz` writes to `kop.d_Zs` (not `D.d_Zs`), so calling it
between `precalc_grad` and `grad` does not corrupt the precomputed $R_m\xi$ values.
`grad_deviceAsyncLaunch` overwrites `D.d_Zs[0]` and `D.d_Ys[0]` as temporaries each
call; `D.d_Zs[m>0]` and `D.d_Ys[m>0]` are only read, never written, per call.

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 5 (DEPRECATED) — Free-field lattice divergence

**ABANDONED (2026-05-28).** Tested with both spatial and temporal links; max|div|=154 on
the free-field background. Root cause not identified; approach superseded by Chunk 5b
(stochastic expectation-value Ward identity).

---

### Chunk 5b — Stochastic Ward identity under expectation value

**Physics.** The exactly conserved Ward identity eq. (3.32) implies

$$\sum_{z:\text{nn of }w} J_{V,+}^{wz} = 0$$

where, for $N_f$ degenerate flavors, the complex-valued current is

$$J_{V,+}^{wz} = -N_f\,\text{tr}\!\left(K^{wz} D_\text{ov}^{-1}\right).$$

Both the real and imaginary parts must sum to zero over all neighbors of $w$.

Fix $w = (s=0, ix=0)$ and sum over all neighbors $z$ (both spatial and temporal links incident to $w$).

**Stochastic estimator.** Use $N_\text{hits}$ Z2$\times$Z2 noise vectors $\eta_k$ (one noise at every lattice site, no dilution):

$$J_{V,+}^{wz} \approx -N_f \cdot \frac{1}{N_\text{hits}} \sum_{k=1}^{N_\text{hits}} \eta_k^\dagger\, K^{wz}\, \phi_k, \quad D_\text{ov}\, \phi_k = \eta_k.$$

The per-hit estimator $\eta_k^\dagger K^{wz}\phi_k$ is complex; multiply by $-N_f$ to get $J_{V,+}^{wz}$. Report Re and Im separately.

Forward solve $D\phi = \eta$ via normal equations: $D^\dagger D\phi = D^\dagger\eta$. No adjoint term or adjoint solve needed.

**Z2$\times$Z2 noise.** Each component: real part $\in \{+1/\sqrt{2}, -1/\sqrt{2}\}$ and imaginary
part $\in \{+1/\sqrt{2}, -1/\sqrt{2}\}$ independently with equal probability. Satisfies
$\mathbb{E}[\eta_i \overline{\eta_j}] = \delta_{ij}$.

**Gauge background for Chunk 5b:** free-field $U=1$ (default-constructed). The Ward identity
holds exactly for any single gauge background (Ginsparg-Wilson guarantee); the stochastic
estimator approximates the all-to-all trace and should converge to zero with $N_\text{hits}$.

**Solver setup** (only DH and DHD wrappers needed):

```cpp
auto f_DH  = std::bind(&Fermion::adj_deviceAsyncLaunch,  &D, _1, _2);
auto f_DHD = std::bind(&Fermion::DHD_deviceAsyncLaunch,  &D, _1, _2);
LinOpWrapper M_DH(f_DH), M_DHD(f_DHD);
MatPoly op_DH;  op_DH.push_back(cplx(1.0), {&M_DH});
MatPoly op_DHD; op_DHD.push_back(cplx(1.0),{&M_DHD});
```

**Forward solve** $D\phi = \eta$ (normal equations):
```cpp
op_DH.from_cpu<N>(DH_eta.field, eta.field);   // DH_eta = D^dag eta
op_DHD.solve<N>(phi.field, DH_eta.field, Comp::TOL_OUTER);
```

**Noise generation** (`includes/valence_claude.h:121`): `eta.fill_z2_source(rng)` — all-to-all
Z2$\times$Z2 at every site and spinor component; satisfies $\mathbb{E}[\eta_i\bar\eta_j]=\delta_{ij}$.

**Analytic reference values** (small lattice only, guarded by `if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4)`):
Loop over all $N$ basis vectors $e_k$; for each solve $D\phi_k = e_k$ and accumulate
$(K^{wz}\phi_k)[k]$ (full complex) per link. After the loop multiply by $-N_f$ to get
$J_{V,+}^{wz}$. Prints before the stochastic hit loop; final summary shows
`Re_analy`, `Re_stoch`, `err_Re`, `sigma_Re` and Im counterparts.

**Procedure.**

1. Use the default-constructed free-field $U$ (already set up in Chunk 2).
2. Set up solver wrappers (above); allocate `FermionVector eta, DH_eta, phi`.
3. (Optional) Analytic pass: basis-vector loop over $k=0\ldots N-1$ to compute `JV_exact[j]` and `JV_exact_t[2]`. Guarded by `if constexpr` for small lattice.
4. For each hit $k = 1,\ldots,N_\text{hits}$:
   a. `eta.fill_z2_source(rng)` — draw new noise.
   b. Forward solve: compute $\phi_k = D^{-1}\eta_k$ via normal equations.
   c. For each neighbor $z$ of $w$ (spatial and temporal links incident to $w$):
      - `kop.apply_k(d_k_phi, d_phi_dev, U, el)` $\to$ $K^{wz}\phi_k$
      - $C_k^{wz} = \text{Im}(\eta_k^\dagger K^{wz}\phi_k)$; store per-hit (real).
      - Temporal backward link: include sign $-1$.
5. After $N_\text{hits}$: $J_{V,+}^{wz} \approx -N_f \cdot \text{mean}(C^{wz})$ (complex).
   Report `Re_analy`, `Re_stoch`, `err_Re`, `sigma_Re`, and Im counterparts per link and for the divergence sum (should be $\approx 0$).

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 6 (FINAL) — Stochastic Ward identity on thermalized configurations (interacting theory)

**Physics.** For the interacting theory the Ward identity is

$$\sum_{z \sim w} \langle J_V^{wz} \rangle = 0,$$

where $\langle \cdot \rangle$ denotes the full path-integral average over both fermions and gauge fields. The Ginsparg-Wilson relation guarantees exact conservation for every gauge background $U$, so the identity holds **config-by-config**:

$$\sum_{z \sim w} \text{tr}\!\left(K^{wz} D_\text{ov}^{-1}[U]\right) = 0 \quad \forall\, U$$

(both Re and Im parts). This is strictly stronger than a statement about the gauge average: the divergence is zero for each $U$ individually, not just on average. Chunk 6 verifies this under realistic gauge fluctuations where $K^{wz}$ and $D_\text{ov}^{-1}$ both depend non-trivially on $U$.

**Stochastic estimator.** Identical to Chunk 5b: $N_\text{hits}$ Z2$\times$Z2 noise vectors, one forward solve $D[U]\phi_k = \eta_k$ per hit, $J_{V,+}^{wz} \approx -N_f \cdot \text{mean}(\eta^\dagger K^{wz}\phi)$ (full complex). The stochastic error shrinks as $1/\sqrt{N_\text{hits}}$; the exact answer is zero.

**Config loading.** Follow `disc_claude.cu:204--257`: discover data directories matching the run parameters, load each binary trajectory file, call `D.update(U)` to refresh the Overlap operator for the new gauge background, then run the estimator.

**Output format:** one line per config:
```
config <k>   div_estimate= <val>   nhits=16
```

**Files:** `check_conserved_current_claude.cu`

---

## 5. Open Questions

1. **`base.links` type** — confirmed as `std::vector<std::array<Idx,2>>` via
   `gauge_ext.h:318`. Iteration `for(const auto& lk : base.links)` gives `BaseLink`
   directly.

## 6. Resolved Implementation Details (Chunk 5b)

- **No adjoint solve needed**: estimator is $\eta^\dagger K^{wz}\phi$ (full complex) with one forward
  solve $D\phi=\eta$; no adjoint solve. Correct formula: $J_{V,+}^{wz} = -N_f\,\text{tr}(K^{wz}D^{-1})$ (complex; both Re and Im sum to zero by Ward identity). Code reports $-N_f \cdot \text{mean}(\eta^\dagger K\phi)$ with Re and Im tracked separately.
- **Noise source**: `FermionVector::fill_z2_source(rng)` in `valence_claude.h:121` — all-to-all
  Z2$\times$Z2 at every site, no dilution.
- **Gauge background**: Chunk 5b uses free-field $U=1$; Chunk 6 loops over thermalized configs.
- **Temporal sign**: spatial links always sign $+1$ (directed $\{ix, nn\}$); temporal forward link
  sign $+1$; temporal backward link $\{(s-1+N_t)\%N_t, ix\}$ sign $-1$.
- **Analytic reference**: basis-vector loop over $N$ vectors; guarded by
  `if constexpr (Comp::N_REFINE == 1 && Comp::Nt <= 4)` for runtime safety.

---

## 6. Resolved Decisions

- $K^{wz}$ is NOT Hermitian; Hermiticity test removed. Proof: $(K^{wz})^\dagger \neq K^{wz}$ from eq. (3.34) — Term A flips $W^{wz}$ to $W^{wz\dagger}$ on the wrong side, Term B moves $X$ from left to $X^\dagger$ on right.
- Ward identity is checked as a vector (full site-by-site divergence), not a scalar.
- Check 2 (Ward identity on configs) is the final and most important check.
- Config loading follows `disc_claude.cu` pattern exactly.
