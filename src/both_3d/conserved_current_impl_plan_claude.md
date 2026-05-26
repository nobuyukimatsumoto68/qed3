# Plan: Exactly Conserved Vector Current for Overlap Fermions (Zolotarev)

## Purpose

Implement the two-component overlap current kernel operator $k^{wz}$ for spatial links,
derived from the Zolotarev rational approximation of the sign function (top-left block of
eq. 3.15 in `qed3int-1.pdf`).  For a given right source $\xi$, the operator returns the
device vector $k^{wz}\xi$; the caller contracts it with stochastic left vectors (one-end
trick in `meson_pq_wall_v2_claude.cu`).  This plan covers only the header
`conserved_current_claude.h`; the measurement driver is out of scope.

---

## Physics Background

### Two-component overlap operator

The code operates entirely in the two-component sector (`NS=2`).  The two-component
overlap operator is (Appendix B.3 of `qed3int_v2-1.pdf`, eq. B.3):

$$D_\text{ov} = 1 + V, \qquad V = \frac{D_W}{\sqrt{D_W^\dagger D_W}},$$

where $D_W$ is the two-component Wilson operator with negative mass $M$ (eq. B.1).
With the Zolotarev rational approximation and partial-fraction decomposition (eq. 3.6 of
`qed3int-1.pdf`):

$$D_\text{ov} = 1 + X\!\left\{p_0 + \sum_{\ell=1}^{\ell_\text{max}} p_\ell (X^\dagger X + q_\ell)^{-1}\right\},$$

where $X \equiv D_W / \lambda_M$ is the scaled two-component Wilson operator.  In the code:
`A[m]` $= p_\ell$, `cp[m]` gives $q_\ell = -k^2/c'_m$, `C` is the $p_0$-related prefactor.

### Vector variation of the action

The two-component action for one flavor pair $(\xi, \eta)$ is (eq. 1.36 of `qed3int_v2-1.pdf`):

$$S = \eta^\dagger D_\text{ov}\,\xi + \xi^\dagger D_\text{ov}^\dagger\,\eta.$$

We focus on the first term $\eta^\dagger D_\text{ov}\xi$.  Under a site-local U(1) rotation with parameter $\theta_w$ (eq. 3.5 of `qed3int-1.pdf`):

$$\delta_V S = \eta^\dagger_x [D_\text{ov}, \Theta]_{xy} \xi_y,$$

where $\Theta_{xy} = \theta_x \delta_{xy}$.  Since only nearest-neighbor links contribute
nontrivially to $[X, \Theta]$ (eq. 3.8):

$$[X, \Theta]_{xy} = \delta\theta_{xy}\, W_{xy}, \qquad \delta\theta_{xy} \equiv \theta_x - \theta_y.$$

Introducing the redundant matrix index notation (eq. 3.12--3.13):

$$\mathcal{W}_{xy}^{wz} \equiv \delta_{xw}\delta_{yz} W_{wz} - \delta_{xz}\delta_{yw} W_{zw},$$

$$[X, \Theta] = \sum_{\langle w,z\rangle} \delta\theta_{wz}\, \mathcal{W}^{wz}.$$

The action variation becomes (eq. 3.14):

$$\delta_V S = \sum_{\langle w,z\rangle} \sum_{x,y} \delta\theta_{wz}\, \eta^\dagger_x (k^{wz})_{xy} \xi_y.$$

### Two-component current kernel $k^{wz}$ (top-left block of eq. 3.15)

$$k^{wz} = \mathcal{W}^{wz}\!\left\{p_0 + \sum_{\ell=1}^{\ell_\text{max}} p_\ell (X^\dagger X + q_\ell)^{-1}\right\}$$
$$- X \sum_{\ell=1}^{\ell_\text{max}} p_\ell (X^\dagger X + q_\ell)^{-1}\!
\left(-\mathcal{W}^{wz\dagger} X + X^\dagger \mathcal{W}^{wz}\right)
(X^\dagger X + q_\ell)^{-1}.$$

Here $X = D_W/\lambda_M$ is the two-component scaled Wilson operator and $\mathcal{W}^{wz}$
is the two-component antisymmetric current kernel (built from $W_{wz}$, see below).

Denote the shifted resolvent for pole $\ell$ as $R_\ell \equiv (X^\dagger X + q_\ell)^{-1}$
and define the resolvent vectors for source $\xi$:

$$Z_\ell \equiv R_\ell \xi, \qquad Z_0 \equiv \xi + \sum_\ell A[\ell]\, Z_\ell.$$

The vector $k^{wz}\xi$ is then:

$$k^{wz}\xi = C\!\left[\mathcal{W}^{wz} Z_0
- X \sum_{\ell=1}^{\ell_\text{max}} A[\ell]\, R_\ell\!\left(X^\dagger \mathcal{W}^{wz} Z_\ell
  - \mathcal{W}^{wz\dagger} X Z_\ell\right)\right].$$

### Two-component Wilson current kernel $W_{wz}$ (eq. B.1--B.2 of `qed3int_v2-1.pdf`)

The code is always two-component (`NS=2`).  The two-component Wilson current kernel on
link $(w,z)$ is obtained from the hopping term of $D_W$ (eq. B.1) by differentiating
with respect to the gauge phase $\theta_{wz}$:

$$W_{wz} \equiv -\kappa_{wz}\, P_{wz}\, U_{wz}, \qquad P_{wz} \equiv \frac{1}{2}(1 - e^a_{wz}\sigma_a)\,\Omega_{wz},$$

where $e^a_{wz}$ is the vielbein along link $(w,z)$, $\sigma_a$ are Pauli matrices, and
$\Omega_{wz}$ is the spin connection matrix.  This is the two-component projector $P_{wz}$
from eq. B.2 of `qed3int_v2-1.pdf`.

$$W_{wz} \equiv -\kappa_{wz}\, P_{wz}\, U_{wz}, \qquad P_{wz} \equiv \frac{1}{2}(1 - e^a_{wz}\sigma_a)\,\Omega_{wz}.$$

In code, `DiracExt::d_coo_format` (`dirac_ext.h:379`) computes the hopping matrix entry
with an extra factor $i$ (from $\partial U/\partial\theta$ in the force):

```cpp
// TEMPLATE: dirac_ext.h:379  (pos branch; includes factor I* for force, not for current)
const MS tmp = 0.5 * bd.kappa[il] * ( -r*sigma[0] + bd.gamma(ix, iy) )
             * I * std::exp( I * u.sp(s, BaseLink{ix,iy}) ) * bd.Omega(ix, iy);
```

For the current kernel we remove the $i$ factor.  The two-component $W_{wz}$ (pos branch,
row $w$, col $z$) and $W_{zw}$ (neg branch, row $z$, col $w$) are therefore:

```cpp
// pos branch: W_{wz}, inserted at row Nx*s+NS*ix, col Nx*s+NS*iy
const MS W_pos = 0.5 * bd.kappa[il] * ( -r*sigma[0] + bd.gamma(ix, iy) )
               * std::exp( I * u.sp(s, BaseLink{ix,iy}) ) * bd.Omega(ix, iy);

// neg branch: -W_{zw}, sign absorbed as per eq. 3.12
// d_coo_format neg has leading minus on MS tmp; same convention kept here
const MS W_neg = -0.5 * bd.kappa[il] * ( -r*sigma[0] + bd.gamma(iy, ix) )
               * std::exp( I * u.sp(s, BaseLink{iy,ix}) ) * bd.Omega(iy, ix);
```

The antisymmetric combination $\mathcal{W}^{wz}$ (eq. 3.12) has nonzero entries at
(row $w$, col $z$) with value $W_{wz}$ and (row $z$, col $w$) with value $W_{zw}$
(see sign structure at `dirac_ext.h:393`):

$$(\mathcal{W}^{wz})_{xy} = \delta_{xw}\delta_{yz} W_{wz} - \delta_{xz}\delta_{yw} W_{zw}.$$

Note: `W_neg` already carries the minus sign from eq. 3.12 via the code convention
at `dirac_ext.h:393`, so the COO entries for $\mathcal{W}^{wz}$ follow the pos branch
directly for the $(w,z)$ block and the neg branch directly for the $(z,w)$ block.

### Exact lattice conservation (eq. 3.18--3.19)

$$\delta_V S = \sum_w \theta_w \sum_{z:\, \text{nn of } w} j^{wz} = 0 \quad \forall \theta_w
\implies \sum_{z:\, \text{nn of } w} j^{wz} = 0.$$

This is the Ward identity used for the numerical check.

---

## References

- N. Matsumoto et al. (internal notes), "qed3int-1.pdf", QED3 project document, 2025 --
  Section 3 ("Exactly conserved lattice currents"), derivation of $\mathcal{K}^{wz}$ and
  the Ward identity check.
- N. Matsumoto, "Notes on QED3 (interacting theory)", "qed3int_v2-1.pdf", 2026-05-25 --
  Section 1 (parity and flavor symmetry, two mass types); Section 2 (pseudofermion
  action eq. 2.3, reweighting factor eq. 2.4): resolves open question \#4 on mass convention.
- N. Matsumoto et al. (internal notes), "qed3int-2.pdf", QED3 project document, 2026 --
  Section 3.1 (four-component derivation, eqs. 3.1-3.19); Section 3.2 (two-component
  block structure, eqs. 3.20-3.22): resolves open question \#1 on two-component form.
- R. Narayanan and H. Neuberger, "Chiral fermions on the lattice", Phys. Lett. B 302 (1993)
  -- original overlap fermion definition.
- A.D. Kennedy, "Algorithms for dynamical fermions", hep-lat/0402038 (2004) -- Zolotarev
  rational approximation and elliptic function coefficients (cited in `overlap_wmass_claude.h:55`).

---

## Symbol Dictionary

| Symbol | Meaning | Code object |
|---|---|---|
| $D_\text{ov}$ | Two-component overlap Dirac operator | `OverlapWMass<WilsonDirac>` (`overlap_wmass_claude.h`) |
| $D_W$ | Two-component Wilson operator with negative mass | `DiracExt<Base, DiracS2Simp>` (`dirac_ext.h`) |
| $X = D_W / \lambda_M$ | Scaled two-component Wilson operator | `(1/lambda_max) * M_DW` (`overlap_wmass_claude.h:330`) |
| $\lambda_M$ | Largest spectral value of $D_W$ | `lambda_max` (member of `OverlapWMass`) |
| $p_0$ | Constant term coefficient (rational approx) | `C` (member of `Zolotarev`) |
| $p_\ell$ | Partial-fraction numerator residues | `A[m]`, $m = 1, \ldots, \ell_\text{max}$ (`Zolotarev`) |
| $q_\ell$ | Partial-fraction pole offsets | `-k*k/cp[m]` (`Zolotarev`) |
| $R_\ell = (X^\dagger X + q_\ell)^{-1}$ | Shifted resolvent, pole $\ell$ | result of `Op.solveAsync<N>` in `mult_deviceAsyncLaunch` |
| $P_{wz} = \frac{1}{2}(1 - e^a_{wz}\sigma_a)\Omega_{wz}$ | Two-component Wilson projector (eq. B.2) | `0.5*(-r*sigma[0] + bd.gamma(ix,iy))*bd.Omega(ix,iy)` (`dirac_simp.h:269`) |
| $\sigma_a$ ($a=1,2,3$) | Pauli matrices (two-component spin basis) | implicit in `bd.gamma(ix,iy)` (`dirac_simp.h`) |
| $W_{wz} = -\kappa_{wz} P_{wz} U_{wz}$ | Two-component Wilson current kernel on link $(w,z)$ | pos branch of `DiracExt::d_coo_format` with factor $I*$ removed (`dirac_ext.h:379`) |
| $\mathcal{W}^{wz}$ | Antisymmetric two-component current operator (eq. 3.12) | new `COO<N>` built per-link (analogous to `COO` in `grad_deviceAsyncLaunch`) |
| $k^{wz}$ | Two-component overlap current kernel (eq. 3.15, top-left block); maps spinor to spinor | `apply_k_wz` in `ConservedCurrent` |
| $Z_\ell = R_\ell \xi$ | Resolvent vector for pole $\ell$; precomputed in `precalc_src` | `D.d_Zs[m]`, $m=1,\ldots,\ell_\text{max}$ |
| $Z_0 = \xi + \sum_\ell A[\ell] Z_\ell$ | Accumulated source vector | `D.d_Zs[0]` |
| $\xi$ | Source spinor (two-component, $2N$ entries) | `CuC* d_xi` on device |
| $\kappa_{wz}$ | Link hopping coefficient | `bd.kappa[il]` (`DiracS2Simp`); `kappa_t[ix]` for temporal links |
| $e^a_{wz}$, $\Omega_{wz}$ | Vielbein direction and spin-connection matrix | `bd.gamma(ix,iy)`, `bd.Omega(ix,iy)` |
| $U_{wz} = e^{i\theta_{wz}}$ | U(1) gauge link phase | `u.sp(s, BaseLink{ix,iy})` (spatial) / `u.tp(s, ix)` (temporal) |

---

## Scope and Constraints

**Target file (new):**
- `src/both_3d/includes/conserved_current_claude.h` — `ConservedCurrent` struct with three functions

**Existing files read (not modified):**
- `src/both_3d/includes/overlap_wmass_claude.h` — `OverlapWMass` (resolvent solve pattern)
- `src/both_3d/includes/dirac_ext.h` — `DiracExt::d_coo_format` (COO build pattern)
- `src/both_3d/includes/sparse_matrix.h` — `COO<N>`, `COOEntry`, `do_it()`, `do_conjugate()`
- `src/both_3d/includes/matpoly.h` — `MatPoly`, `solveAsync`, `on_gpuAsync`, `dotAsync`
- `src/both_3d/includes/dirac_simp.h` — `DiracS2Simp` (kappa, gamma, Omega)

**Out of scope:**
- Temporal links (spatial only)
- Any four-component objects: the code is always two-component (`NS=2`)
- Measurement driver, output, correlators, Ward-identity checks
- New CUDA kernels (reuse `Taxpy_gen`, `MatPoly`, `COO`)

---

## Data Structures

### `ConservedCurrent` struct (host-resident)

```
struct ConservedCurrent<OverlapOp>
  const OverlapOp& D           // const ref to OverlapWMass; read-only access to
                               //   M_DW, M_DWH, A[], cp[], k, C, lambda_max, size,
                               //   stream[], handle[], nstreams
  std::vector<CuC*> d_Zs      // d_Zs[0]=Z_0; d_Zs[1..size-1]=Z_ell; owned here
  CuC* d_tmp1, d_tmp2, d_tmp3 // per-link scratch (size N each) for apply_k_wz
```

Device memory (`d_Zs`, scratch) allocated in constructor, freed in destructor,
mirroring the pattern in `overlap_wmass_claude.h:167-196`.

---

## Implementation Chunks

### Chunk 1 — `build_W_wz(COO<N>& coo, int s, BaseLink lk)`

Build the antisymmetric two-component COO sparse matrix $\mathcal{W}^{wz}$ for one
spatial link.  Entries:

$$(\mathcal{W}^{wz})_{xy} = \delta_{xw}\delta_{yz} W_{wz} - \delta_{xz}\delta_{yw} W_{zw}$$

with the two-component Wilson current kernel:

$$W_{wz} = -\kappa_{wz}\,P_{wz}\,U_{wz}, \qquad P_{wz} = \tfrac{1}{2}(1 - e^a_{wz}\sigma_a)\,\Omega_{wz}.$$

Template: `dirac_ext.h` pos/neg branches of `d_coo_format`, with the factor `I*` removed.
Row/col indexing: `Nx*s + NS*ix + {0,1}`.  Call `coo.do_it()` after filling entries.
Public helper; also called internally by `apply_k_wz`.

### Chunk 2 — `template<typename Gauge> void apply_k_wz(CuC* d_result, const CuC* d_xi, const Gauge& U, int s, BaseLink lk)`

Compute the resolvent vectors and the current vector in one call.
Builds both the forward and adjoint COOs internally from `(U, s, lk)`; caller does not
manage COO objects.

**Step 1 — Z solves** (parallel CG, mirrors `precalc_grad_deviceAsyncLaunch:409-418`):

```cpp
#pragma omp parallel for num_threads(D.nstreams)
for(int m=1; m<D.size; m++) {
    const int istream = omp_get_thread_num();
    MatPoly Op(D.handle[istream], D.stream[istream], istream);
    Op.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
    Op.push_back(cplx(-D.k*D.k/D.cp[m]), {});
    Op.solveAsync<N>(d_Zs[m], d_xi, Comp::TOL_INNER);
    CUDA_CHECK(cudaStreamSynchronize(D.stream[istream]));
}
CUDA_CHECK(cudaDeviceSynchronize());
// Z_0 = d_xi + sum_m A[m] Z_m
CUDA_CHECK(cudaMemcpy(d_Zs[0], d_xi, N*CD, D2D));
for(int m=1; m<D.size; m++)
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_Zs[0], D.A[m], d_Zs[m], d_Zs[0]);
CUDA_CHECK(cudaDeviceSynchronize());
```

**Step 2 — build COOs**: call `build_W_wz(coo_W, U, s, lk)` for the forward matrix.
For the adjoint `coo_WH`, call `D.DW.d_coo_format` again and manually conjugate entries
and swap `(i,j)`, then call `coo_WH.do_it()`.  Reason: `COO::do_conjugate()`
(`sparse_matrix.h:192`) has no bounds check on the row loop and would crash on the sparse
$\mathcal{W}^{wz}$ whose adjoint has mostly-empty rows.

Adjoint entry transform: $W^\dagger_{ji} = \overline{W_{ij}} = \overline{-I z} = I\bar{z}$,
where $z$ is the entry from `d_coo_format`.  In code: `e.v = cplx(Complex(imag(z), real(z)))`,
`swap(e.i, e.j)`.  Both `coo_W` and `coo_WH` are local variables; GPU memory is freed
at end of `apply_k_wz`.

**Step 3 — apply kernel**: compute

$$k^{wz}\xi = C\!\left[\mathcal{W}^{wz} Z_0 - X\sum_{\ell=1}^{\ell_\text{max}} A[\ell]\, R_\ell\!\left(X^\dagger \mathcal{W}^{wz} Z_\ell - \mathcal{W}^{wz\dagger} X Z_\ell\right)\right]$$

**Term A** (no CG): `coo_W(d_tmp1, d_Zs[0])`; `d_result = C * d_tmp1`.

**Term B** per pole $\ell$ (one synchronous CG solve per pole):

```
coo_W(d_tmp1, d_Zs[m])                  // d_tmp1 = W Z_ell
XH.on_gpu<N>(d_tmp2, d_tmp1)            // d_tmp2 = X^dag W Z_ell
X.on_gpu<N>(d_tmp3, d_Zs[m])            // d_tmp3 = X Z_ell
coo_WH(d_tmp1, d_tmp3)                  // d_tmp1 = W^dag X Z_ell
d_tmp2 -= d_tmp1                         // d_tmp2 = u_ell
Op.solve<N>(d_tmp1, d_tmp2, TOL_INNER)  // d_tmp1 = R_ell u_ell
X.on_gpu<N>(d_tmp3, d_tmp1)             // d_tmp3 = X R_ell u_ell
d_result -= C * A[m] * d_tmp3           // accumulate
```

Poles processed sequentially (shared scratch `d_tmp1/2/3`).

Returns nothing; `d_result` holds $k^{wz}\xi$.  Intended use: one-end trick stochastic
estimator in `meson_pq_wall_v2_claude.cu` and `disc_claude.cu`, replacing `mult_sigma`.

---

## Resolved Decisions

- **Two-component only**: $k^{wz}$ is the top-left block of eq. 3.15; no four-component objects.
- **Spatial links only**: temporal link overload deferred.
- **Mass prefactor**: $k^{wz}$ itself is massless; caller applies $(1+m^*)$ as needed.
- **Vector output**: `apply_k_wz` returns a device vector $k^{wz}\xi$, not a scalar bilinear.
  The caller (one-end trick in `meson_pq_wall_v2_claude.cu`) contracts with stochastic left vectors.
- **No precalc function**: Z solves and kernel application are combined in `apply_k_wz`;
  `d_Zs[]` is pre-allocated in the struct to avoid `cudaMalloc`/`cudaFree` in the hot path.
- **Per-link CG**: Term B of `apply_k_wz` requires one CG solve per Zolotarev pole per link call.
  Poles are processed sequentially (shared scratch buffers `d_tmp1/2/3`).
- **Adjoint COO**: `COO::do_conjugate()` (`sparse_matrix.h:213`) has no bounds check in its
  row-pointer loop and crashes when the adjoint matrix has empty trailing rows.  Instead the
  adjoint is built by re-calling `d_coo_format`, manually transforming entries
  (`e.v = cplx(Complex(imag(z),real(z)))`, `swap(e.i,e.j)`), and calling `do_it()`.
