# Plan: Implicit Restarted Lanczos + Chebyshev Filtering for QED3 Eigenvalue Computation

## Physics / Goal Summary

Compute the $k = 50$ smallest eigenvalues $|\lambda_i|^2$ of the normal operator
$$A \equiv (D_{\text{ov}} + m)^\dagger (D_{\text{ov}} + m),$$
where $D_{\text{ov}}$ is the Zolotarev overlap Dirac operator with complex mass $m$.
$A$ is Hermitian positive semi-definite; its spectrum is real and non-negative.

The existing `eig_wmass_claude.cu` builds the full dense matrix and calls `cusolverDnXgeev`,
which costs $O(N^3)$ time and $O(N^2)$ memory.  For $N = 3072$ ($L=1$) this is borderline;
for $N = 10752$ ($L=2$) it is infeasible.

The new program uses:
1. **Chebyshev filtering** -- replaces $A$ by $p(A) \equiv T_d(q(A; \alpha, \beta))$ where
   $$q(\lambda; \alpha, \beta) \equiv \frac{2\lambda^2 - (\alpha^2 + \beta^2)}{\alpha^2 - \beta^2}, \quad \alpha, \beta > 0.$$
   The polynomial $T_d$ (Chebyshev of first kind, degree $d$) maps $[\beta, \alpha]$ to $[-1, 1]$
   and amplifies eigenvalues with $|\lambda| < \beta$ relative to those with $|\lambda| > \alpha$.
   Parameters: `alpha` (upper edge of target window) and `beta` (spectral upper bound) are
   positional runtime arguments.
2. **Implicitly Restarted Lanczos (IRL)** -- builds a Krylov space of dimension $m = 150$,
   extracts $k = 50$ wanted eigenpairs, and deflates using shifted-QR on the $m \times m$
   tridiagonal matrix.  Each restart reuses the $k$-dimensional compressed Lanczos relation
   (equations 2.17 -- 2.26 in `lanczos.pdf`).

Output: $k = 50$ eigenpairs printed to `stdout`/`clog` in the same format as `eig_wmass_claude.cu`:
```
// one line per eigenpair: i  real(lambda)  imag(lambda)  abs(lambda)
```
(eigenvalues of $A$ are real, so `imag` will be 0.)

## All Files to Be Modified or Created

| File | Role |
|------|------|
| `includes/lanczos_claude.h` (NEW) | Self-contained IRL + Chebyshev header; no Grid dependency |
| `eig_lanczos_claude.cu` (NEW) | Entry point: `main()`; mirrors `eig_wmass_claude.cu` structure |

No existing files are modified.

## References

- N. Matsumoto, "Large-scale eigenvalue problem" (notes, May 2026) -- `/mnt/barracuda22/qed3/lanczos.pdf` -- IRL algorithm (Sec. 2), Chebyshev filtering (Sec. 3), notation used throughout this plan.
- Z. Bai et al., "Templates for the solution of algebraic eigenvalue problems", SIAM 2000 -- underlying IRL theory.
- P. Boyle et al., Grid library `ImplicitlyRestartedLanczos.h`, `/mnt/baracuda_14/grid_claude/Grid/Grid/algorithms/iterative/ImplicitlyRestartedLanczos.h` -- reference implementation consulted for `step()`, `QR_decomp()`, `diagonalize_Eigen()`, and convergence test.

## Resolved Questions

| # | Question | Answer |
|---|----------|--------|
| Q1 | Which operator? | $(D_{\text{ov}}+m)^\dagger (D_{\text{ov}}+m)$ via `DHD_deviceAsyncLaunch`; wrapped by `LinOpDHDWrapper` |
| Q2 | How many eigenpairs / Krylov dim? | $k = 50$ converged, $m = 150$ Krylov vectors |
| Q3 | Output format? | Same as `eig_wmass_claude.cu`: one line per pair, `i real imag abs` to `stdout`/`clog` |
| Q4 | Chebyshev range? | `alpha` and `beta` as positional argv[3] and argv[4] |

## Ordered Implementation Chunks

---

### Chunk 1 -- `eig_lanczos_claude.cu`: boilerplate and `Comp` namespace

**Description.** Copy the header block of `eig_wmass_claude.cu` verbatim (includes, type
aliases, `Comp` namespace, `dir` string, GPU setup).  Add argv parsing for `mass_re`,
`mass_im`, `alpha`, `beta`.  Construct `Base`, `WilsonDirac`, `Gauge`, `Fermion`
(`OverlapWMass<WilsonDirac>`), call `Dov.update(U)` on a cold (unit) gauge.  Instantiate
`LinOpDHDWrapper<Fermion> A_op(Dov)`.

**Files:** `eig_lanczos_claude.cu`

**Template** (`eig_wmass_claude.cu:32--186`):
```cpp
// TEMPLATE
namespace Comp{
  constexpr int N_REFINE = 4;
  constexpr int Nt = 24;
  ...
}
...
Fermion Dov(DW, mass, 21);
Dov.update(U);
LinOpDHDWrapper<Fermion> A_op(Dov);  // applies (D+m)^dag (D+m)
```

**Notes.**
- `Comp::N_REFINE` and `Comp::Nt` are the only compile-time knobs the user will change between lattice sizes; keep them identical to `eig_wmass_claude.cu`.
- argv order: `[mass_re] [mass_im] [alpha] [beta] [cheb_degree]` (degree defaults to 12 if omitted).

---

### Chunk 2 -- `includes/lanczos_claude.h`: Chebyshev filter as a `LinOp`

**Description.** Define `struct ChebyshevFilter : public LinOp`.  Given the base operator
$A$ (`LinOp&`), degree $d$, and parameters `alpha`, `beta`, it applies

$$p(A) v = T_d(q(A;\alpha,\beta)) v$$

via the three-term Chebyshev recurrence (no explicit polynomial coefficients needed):

$$v_0 = v, \quad v_1 = q(A;\alpha,\beta) v,$$
$$v_{j+1} = 2\, q(A;\alpha,\beta) v_j - v_{j-1}, \quad j = 1, \ldots, d-1,$$

where applying $q(A;\alpha,\beta)$ to a vector $w$ means:
$$q(A;\alpha,\beta) w = \frac{2 A w - (\alpha^2 + \beta^2) w}{\alpha^2 - \beta^2}.$$

Each application of $q(A;\alpha,\beta)$ costs exactly one application of the base operator $A$
plus two cuBLAS `Zdaxpy` calls and one `Zdscal`, all on device pointers.

**Files:** `includes/lanczos_claude.h`

**Template** (`includes/sparse_matrix.h:4--9` for `LinOp` interface):
```cpp
// TEMPLATE
struct LinOp{
  using T = CuC;
  virtual void operator()( T* d_res, const T* d_v ) const = 0;
  virtual void Async( T* d_res, const T* d_v,
                      const cudaStream_t stream) const = 0;
};
```

**Template** (`includes/gpu_header.h:127` for in-place `Taxpy`):
```cpp
// TEMPLATE
void Taxpy(T* d_res, const T d_a, const T* d_x, const T* d_y){
```

**Design decisions.**
- `ChebyshevFilter` owns two temporary `CuC*` device buffers (`d_v0`, `d_v1`) of size
  `Comp::N`, allocated in constructor, freed in destructor.
- The `operator()` entry point calls `d_v0 = v`, then iterates the recurrence in-place.
- No cuBLAS handle stored internally: use a local `MatPoly` (default constructor) for
  `Zdscal` and `dot`; use `Taxpy_gen` kernels directly for the `daxpy` steps, exactly as
  done in `matpoly.h:119`.
- `alpha`, `beta`, `degree` stored as `const double` / `const int` members.

---

### Chunk 3 -- `includes/lanczos_claude.h`: data structures for IRL

**Description.** Define the IRL workspace struct `IRLState` holding all Lanczos vectors and
the tridiagonal representation on the host.

```cpp
// signature (in lanczos_claude.h)
struct IRLState {
  const int Nk;  // wanted eigenpairs (50)
  const int Nm;  // Krylov dimension  (150)
  int  Nconv;    // converged count after calc()

  // device: Nm Lanczos basis vectors, each of length Comp::N
  std::vector<CuC*> d_evec;  // size Nm, each cudaMalloc'd Comp::N*CD

  // host: tridiagonal T_m
  std::vector<double> alpha;  // diagonal,    length Nm
  std::vector<double> beta;   // sub-diagonal, length Nm (beta[Nm-1] unused)

  // host: converged eigenvalues
  std::vector<double> eval;   // length Nm (trimmed to Nk after convergence)

  // host: Eigen Nm x Nm rotation matrix Q (for implicit shifts)
  Eigen::MatrixXd Qt;

  // device: residual vector f (current Lanczos direction)
  CuC* d_f;

  IRLState(int Nk_, int Nm_);
  ~IRLState();
};
```

**Files:** `includes/lanczos_claude.h`

**Template** (`includes/overlap_wmass_claude.h:151--189` for CUDA alloc/dealloc pattern):
```cpp
// TEMPLATE
for(int m=0; m<size; m++) {
  CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
  ...
}
```

**Design decisions.**
- All $m = 150$ Lanczos vectors reside on device (GPU) as `CuC*`; size per vector is
  `Comp::N * sizeof(cuDoubleComplex)`.
- Memory cost: $150 \times 3072 \times 16 \approx 7\,\text{MB}$ for $L=1$;
  $150 \times 10752 \times 16 \approx 26\,\text{MB}$ for $L=2$. Both fit comfortably.
- `alpha`, `beta`, `eval`, `Qt` are CPU-side (Eigen / `std::vector<double>`), consistent
  with the Grid IRL reference which does all tridiagonal arithmetic on CPU.

---

### Chunk 4 -- `includes/lanczos_claude.h`: Lanczos `step()` function

**Description.** Implement one step of the basic Lanczos iteration (equations 2.8 -- 2.13
of `lanczos.pdf`; also Saad p.195 comment in Grid IRL reference, line 449):

Given basis vector `d_evec[k]` and previous beta `beta[k-1]`, compute:
1. `d_f = PolyOp( d_evec[k] )` -- apply filtered operator.
2. If `k > 0`: `d_f -= beta[k-1] * d_evec[k-1]`.
3. `alpha[k] = Re( dot(d_evec[k], d_f) )`.
4. `d_f -= alpha[k] * d_evec[k]`.
5. `beta[k] = sqrt( dot2self(d_f) )`.
6. `d_evec[k+1] = d_f / beta[k]` (if `k < Nm-1`).
7. Periodic full re-orthogonalization against all previous vectors when `k % orth_period == 0` (default `orth_period = 4`).

```cpp
// signature
void lanczos_step(const LinOp& PolyOp,
                  IRLState& s,
                  int k,
                  cublasHandle_t handle);
```

**Files:** `includes/lanczos_claude.h`

**Template** (`includes/matpoly.h:152--175` for `Zdscal` and `dot`):
```cpp
// TEMPLATE
CUBLAS_CHECK( cublasZdscal(handle, N, &alpha, x, 1) );
CUBLAS_CHECK( cublasZdotc(handle, N, x, 1, y, 1, result) );
```

**Template** (`includes/gpu_header.h:127` for `Taxpy`):
```cpp
// TEMPLATE
Taxpy<CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_res, al, d_p, d_x);
```

**Design decisions.**
- All dot products / norms are done via `cublasZdotc`; the real part is taken on host.
- Re-orthogonalization uses sequential `cublasZdotc` then `Taxpy_gen` for each prior vector,
  exactly as `basisOrthogonalize` does in Grid (no block GEMM needed at $m = 150$).
- One shared `cublasHandle_t` is passed in (created in `main()`); no extra cuBLAS handles.

---

### Chunk 5 -- `includes/lanczos_claude.h`: tridiagonal diagonalization and shifted QR

**Description.** Implement two CPU-side routines operating on `std::vector<double>` arrays
and `Eigen::MatrixXd Qt`:

**5a. `diagonalize_Eigen()`** -- uses `Eigen::SelfAdjointEigenSolver` on the $k \times k$
leading block of $T_m$, accumulates eigenvectors into `Qt`. This is a direct port of the
Grid routine at `ImplicitlyRestartedLanczos.h:501--522`.

**5b. `QR_decomp()`** -- one Givens-rotation pass for a single shift $\mu$, updating
`alpha`, `beta`, and `Qt` in place. Direct port of Grid `QR_decomp()` at lines 527 -- 590.

```cpp
// signatures
void diagonalize_Eigen(std::vector<double>& lmd,
                       std::vector<double>& lme,
                       int Nk, int Nm,
                       Eigen::MatrixXd& Qt);

void QR_decomp(std::vector<double>& lmd,
               std::vector<double>& lme,
               int Nk, int Nm,
               Eigen::MatrixXd& Qt,
               double shift, int kmin, int kmax);
```

**Files:** `includes/lanczos_claude.h`

**Template** (`ImplicitlyRestartedLanczos.h:501--590`):
```cpp
// TEMPLATE
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(TriDiag);
for (int i = 0; i < Nk; i++) lmd[Nk-1-i] = eigensolver.eigenvalues()(i);
```

**Design decisions.**
- All tridiagonal operations are CPU-only; the matrix is $150 \times 150$ at most.
- The eigenvalues are sorted in descending order after diagonalization, matching Grid's
  convention (`std::partial_sort(..., std::greater<double>()`); the $p = m - k = 100$
  largest become the implicit shifts.

---

### Chunk 6 -- `includes/lanczos_claude.h`: basis rotation on GPU

**Description.** After the shifted-QR passes accumulate the $m \times m$ rotation `Qt`,
apply `V_k := V_m Q(:, :k)` on the GPU: for each of the $k$ output vectors, compute a
linear combination of the $m$ existing Lanczos vectors using the column of `Qt`.

```cpp
// signature
void basisRotate(std::vector<CuC*>& d_evec,
                 const Eigen::MatrixXd& Qt,
                 int k_lo, int k_hi,   // output range [k_lo, k_hi)
                 int Nm,
                 CuC* d_tmp,           // scratch, size Comp::N
                 cublasHandle_t handle);
```

**Files:** `includes/lanczos_claude.h`

**Template** (`includes/matpoly.h:106--126` for GPU linear combination pattern):
```cpp
// TEMPLATE
for(int i=0; i<vec_mats.size(); i++){
  CUDA_CHECK(cudaMemcpy(d_tmp, d_v0, N*CD, D2D));
  ...
  Taxpy<CuC, N><<<NBlocks, NThreadsPerBlock>>>(d_v, coeffs[i], d_Mv0, d_v);
}
```

**Design decisions.**
- Loop over $j \in [k\_lo, k\_hi)$: zero `d_tmp`, then accumulate
  `d_tmp += Qt(i,j) * d_evec[i]` for $i = 0, \ldots, m-1$ using `Taxpy_gen` with
  `double` coefficient, then `cudaMemcpy(d_evec[j], d_tmp, N*CD, D2D)`.
- After rotation, update the residual vector `d_f` with the scalar from `Qt(Nm-1, Nk-1)`
  and the off-diagonal `beta[Nk-1]`, as in equations 2.21 -- 2.23 of `lanczos.pdf`.
- Scratch buffer `d_tmp` is a single `CuC*` passed in from `main()`.

---

### Chunk 7 -- `includes/lanczos_claude.h`: convergence test and `calc()` driver

**Description.** Implement the top-level `calc()` function that orchestrates the full IRL loop.

```cpp
// signature
void calc(const LinOp& PolyOp,
          const LinOp& HermOp,   // raw A (no filter) for residual check
          IRLState& s,
          const CuC* d_src,      // initial vector (device)
          const double eresid,   // convergence threshold
          const int MaxIter,
          cublasHandle_t handle);
```

Algorithm (following `lanczos.pdf` Sec. 2.3 and Grid `calc()` structure):

1. Normalize `d_src` into `d_evec[0]`.
2. Run $k$ initial Lanczos steps using `PolyOp`.
3. Outer restart loop up to `MaxIter`:
   a. Extend Krylov space from $k$ to $m$ steps.
   b. `diagonalize_Eigen()` on $T_m$; sort eigenvalues; collect $p = m-k$ unwanted shifts.
   c. Apply $p$ shifted-QR passes via `QR_decomp()`.
   d. `basisRotate()` to compress to $k$ vectors; update `d_f` and `beta[k-1]`.
   e. Resume Lanczos from step $k$.
   f. Convergence check: for each of the $k$ candidates compute
      $\rho_j = \Vert A \tilde{u}_j - \tilde{\lambda}_j \tilde{u}_j \Vert / \lambda_{\max}$
      using `HermOp`; accept if $\rho_j < \text{eresid}$.
   g. If all $k$ accepted, exit.
4. On convergence: rotate `d_evec[0..k-1]` by the final `Qt` block; reconstruct exact
   eigenvalues via Rayleigh quotient with `HermOp`.

**Files:** `includes/lanczos_claude.h`

**Template** (`ImplicitlyRestartedLanczos.h:211--445` for overall `calc()` structure).

**Template** (`ImplicitlyRestartedLanczos.h:66--93` for convergence test logic):
```cpp
// TEMPLATE
_HermOp(B, v);
RealD vnum = real(innerProduct(B, v));
RealD vden = norm2(B);
eval = vnum/vden;
v -= eval*B;
RealD vv = norm2(v) / pow(evalMaxApprox, 2.0);
int conv = (vv < eresid*eresid) ? 1 : 0;
```

**Design decisions.**
- `PolyOp` is the `ChebyshevFilter`; `HermOp` is `LinOpDHDWrapper` (raw $A$). They are
  separate to allow the residual check to use the unfiltered operator, exactly as Grid does.
- `evalMaxApprox` is estimated once at startup via 50 power-iteration steps with `HermOp`,
  matching Grid lines 238 -- 258.
- Progress is printed to `std::clog` every restart iteration: iteration number, `beta_k`,
  and the current $k$ Ritz values.

---

### Chunk 8 -- `eig_lanczos_claude.cu`: `main()` wiring and output

**Description.** Wire everything together in `main()`:

1. Parse argv: `mass_re`, `mass_im`, `alpha`, `beta`, `cheb_degree` (default 12).
2. Construct lattice + `Dov` exactly as in `eig_wmass_claude.cu:139--184`.
3. Create `LinOpDHDWrapper<Fermion> HermOp(Dov)`.
4. Create `ChebyshevFilter ChebOp(HermOp, cheb_degree, alpha, beta)`.
5. Allocate `IRLState s(50, 150)`.
6. Build random initial vector `d_src` (use `rng` or a fixed seed).
7. Call `calc(ChebOp, HermOp, s, d_src, 1.0e-8, 200, cublas_handle)`.
8. Print results to `clog` in the format:
   ```cpp
   // for i in 0..Nconv-1:
   std::clog << i << " " << s.eval[i] << " " << 0.0 << " " << s.eval[i] << "\n";
   ```
   (real eigenvalues of $A = (D+m)^\dagger (D+m)$; imag part identically zero.)
9. Free all resources, `cudaDeviceReset()`.

**Files:** `eig_lanczos_claude.cu`

**Template** (`eig_wmass_claude.cu:310`):
```cpp
// TEMPLATE
for(int i=0; i<n; i++)
  std::clog << i << " " << real(W[i]) << " " << imag(W[i]) << " " << abs(W[i]) << std::endl;
```

---

## Open Questions

None that are blocking.  The following are implementation-time judgments left to the coder:

- **Chebyshev degree default**: 12 is a reasonable starting point. Can be changed via argv.
- **`orth_period`**: default 4 (re-orthogonalize every 4 steps). Can be hardcoded.
- **`eresid` default**: `1.0e-8`, matching `Comp::TOL_OUTER`. Can be made an argv parameter in a later iteration.

## Synchronization Convention

After every `Taxpy` or `Taxpy_gen` kernel launch, call `CUDA_CHECK(cudaDeviceSynchronize())` unconditionally.
This applies everywhere in `lanczos_claude.h` and `eig_lanczos_claude.cu` where a custom CUDA kernel is launched.

## Efficiency Decisions Recorded

- Lanczos vectors live entirely on device; only $m \times m$ tridiagonal data crosses to host.
- Single cuBLAS handle shared across all dot-product and scale calls (no per-stream handles in the IRL core).
- Basis rotation uses a serial loop over $k = 50$ output vectors with `Taxpy_gen`; no cublas GEMM needed at this scale.
- Chebyshev recurrence allocates exactly 2 temporary device vectors (persistent in `ChebyshevFilter`), not per-call.
- `MatPoly` is not used inside `ChebyshevFilter`; direct `cublasZdscal` + `Taxpy_gen` calls are used to avoid the `MatPoly` internal `cudaMalloc`/`cudaFree` overhead per application.
