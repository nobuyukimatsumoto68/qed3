# Implementation Plan: Check Program for the Exactly Conserved Vector Current $k^{wz}$

_Created: 2026-05-25_

---

## 1. Physics Motivation

The kernel `apply_k_wz` in `conserved_current_claude.h` implements $k^{wz}$, the
two-component overlap-fermion exactly conserved vector current for a single spatial
link $(s \to \text{lk})$. Its correctness has three independent signatures:

1. **Hermiticity.** The kernel $k^{wz}$ (eq. 3.15 top-left block of `qed3int-1.pdf`)
   must satisfy $\langle \eta, k^{wz} \xi \rangle = \overline{\langle \eta, k^{wz\dagger}\xi \rangle} = \langle k^{wz} \eta, \xi \rangle^*$,
   i.e. $k^{wz}$ is Hermitian. This holds per link and for any gauge background.

2. **Free-field lattice conservation.** For the free-field background ($U = 1$,
   all link phases zero), the summed lattice divergence
   $$(\text{div}\, k\, \xi)(x) = \sum_{\langle w,z\rangle \ni x} \text{sign}(x,wz)\, [k^{wz}\xi](x)$$
   must vanish site-by-site to machine precision for any $\xi$.

3. **Ward identity on thermalized configurations (main check).** For any $\xi$ and
   any thermalized gauge background, the same lattice divergence must vanish to
   $O(\text{TOL\_INNER})$. This is the exact on-lattice Ward identity guaranteed by
   the Ginsparg-Wilson relation.

The COO adjoint construction in `conserved_current_claude.h:86--94` — manual
conjugate-transpose to avoid the bounds-checking bug in `COO::do_conjugate()` — and the
entry transform `e.v = cplx(Complex(imag(z), real(z)))` are the highest-risk code paths
and the primary targets of Check 1.

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
// h_k_xi: host copy of k^{wz}\xi for this link, length Comp::N
// sign convention: site lk[0] rows get +1, site lk[1] rows get -1
// (follows antisymmetry of \mathcal{W}^{wz} and dirac_simp.h:294--338)
void accumulate_divergence(std::vector<Complex>& h_div,
                           const std::vector<Complex>& h_k_xi,
                           int s_lat, const BaseLink& lk);
```

The sign in `accumulate_divergence` must be verified against the COO entries produced
by `build_W_wz` before finalizing (see Open Question 1).

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 4 — Check 1: Hermiticity test

Fill `d_xi` and `d_eta` from host Gaussian random vectors. Then for each `s_lat` in
$[0, N_t)$ and each `lk` in `base.links`:

```
// a = <eta, k^{wz} xi>
kop.apply_k_wz(d_k_xi, d_xi, U, s_lat, lk);
device_to_host_vec(h_k_xi, d_k_xi, N);
a = inner_product(h_eta, h_k_xi);

// b = <k^{wz} eta, xi> = conj(<xi, k^{wz} eta>)
kop.apply_k_wz(d_k_eta, d_eta, U, s_lat, lk);
device_to_host_vec(h_k_eta, d_k_eta, N);
b = conj(inner_product(h_xi, h_k_eta));

// hermiticity residual
print |a - b| / |a|
assert |a - b| / |a| < 1.0e-10
```

No separate adjoint API is needed: the Hermiticity test requires only `apply_k_wz`
called with two different source vectors.

Run on the default-constructed free-field `U`. Optionally repeat on one loaded config.

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 5 — Check 2: Free-field lattice divergence

With `U` at zero-phase (default-constructed):

1. Draw random $\xi$ on host; copy to `d_xi`.
2. Initialize `h_div(N, Complex(0.0, 0.0))`.
3. For `s_lat` in $[0, N_t)$, for `lk` in `base.links`:
   - `kop.apply_k_wz(d_k_xi, d_xi, U, s_lat, lk)`;
   - `CUDA_CHECK(cudaDeviceSynchronize())`;
   - `device_to_host_vec(h_k_xi, d_k_xi, N)`;
   - `accumulate_divergence(h_div, h_k_xi, s_lat, lk)`.
4. Compute `max_i |h_div[i]|` and `norm(h_div)`.
5. Print both; assert `max_i |h_div[i]| < 1.0e-10`.

On the free-field background the CG systems for $Z_\ell$ may still accumulate at
`TOL_INNER`; see Open Question 4 for threshold adjustment.

**Files:** `check_conserved_current_claude.cu`

---

### Chunk 6 (FINAL) — Check 3: Ward identity on loaded configurations

This is the primary and most important check. It validates $k^{wz}$ against
thermalized, interacting gauge configurations read from disk.

**Config discovery:** Follow `disc_claude.cu:204--257`. Build `dir3` from `gsq`, `at`,
`Nf`, `nu0`, `Nt`, `N_REFINE`. Scan `dir3 + "ckpoint_lat." + std::to_string(k)` for
$k = k_\text{ckpt},\, 2k_\text{ckpt},\, \ldots$ using `std::filesystem::exists`.
Process up to 10 configs.

**Per-config loop:**

```cpp
U.read(dir3 + "ckpoint_lat." + std::to_string(k));
D.update(U);
// draw fresh random xi for each config
// h_div = 0
// for s_lat in [0, Nt), for lk in base.links:
//   kop.apply_k_wz(d_k_xi, d_xi, U, s_lat, lk);
//   accumulate_divergence(h_div, ...);
// report max_i |h_div[i]| and norm(h_div) / norm(xi)
// assert max_i |h_div[i]| / norm(xi) < 1.0e-4
```

**Output format:** one line per config:
```
config <k>   max_div= <val>   norm_div/norm_xi= <val>
```

**Threshold rationale:** With `TOL_INNER = 1.0e-9`, $\ell_\text{max} \sim 10$ poles,
$N_\text{links} = 30$, $N_t = 128$: worst-case accumulated CG error
$\sim N_\text{links} \cdot N_t \cdot \ell_\text{max} \cdot \text{TOL\_INNER} \approx 4 \times 10^{-5}$.
Use `1.0e-4 * norm(xi)` as a conservative per-site threshold; report the actual value.

**Files:** `check_conserved_current_claude.cu`

---

## 5. Open Questions

1. **Sign convention in `accumulate_divergence`.** The plan assigns $+1$ to rows
   belonging to site `lk[0]` and $-1$ to rows belonging to site `lk[1]`, based on
   the antisymmetry in `dirac_simp.h:294--338`. However, `d_coo_format` inserts a
   leading minus on the neg branch (`dirac_simp.h:327`). Verify by printing the
   nonzero COO entries of `build_W_wz` for a known free-field link before finalizing
   `accumulate_divergence`.

2. **`base.links` canonical orientation.** `GaugeExt::sp(s, BaseLink{ix,iy})` applies
   `map2sign` ($\pm 1$). When iterating `for(const auto& lk : base.links)`, links are
   already canonical. Confirm that passing `lk` directly to `apply_k_wz` is correct
   with no sign flip needed at the call site.

3. **`d_MemorySets` requirement.** `disc_claude.cu:175` calls
   `d_MemorySets[i].allocate()` before constructing `WilsonDirac`. Confirm whether
   this global is needed in the check program (it is declared in `sparse_dirac.h`).

4. **Free-field divergence tolerance.** On the free-field background the $Z_\ell$ CG
   solves may accumulate error at `TOL_INNER`. If the observed max divergence on
   free-field exceeds $10^{-10}$ (due to CG, not a bug), adjust the Check 2 threshold
   to match and document the finding.

5. **`base.links` type.** Confirmed as `std::vector<std::array<Idx,2>>` via
   `gauge_ext.h:318`. Iteration `for(const auto& lk : base.links)` gives `BaseLink`
   directly.

---

## 6. Resolved Decisions

- Hermiticity test uses two `apply_k_wz` calls; no separate adjoint function needed.
- Ward identity is checked as a vector (full site-by-site divergence), not a scalar.
- Check 3 is the final and most important check; must be the last one in the program.
- Config loading follows `disc_claude.cu` pattern exactly.
