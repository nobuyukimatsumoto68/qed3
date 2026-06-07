# ConservedCurrent -> LinOp refactor (impl plan)

> NOTE (2026-06-04): the operator-wiring follow-up is RESOLVED and lives in the canonical plan
> `conserved_current_correlators_impl_plan_v3_claude.md`, not here. Resolution: the ($++$)
> connected correlator uses $D_m$ on BOTH legs (no massless $D$); only the mechanical cleanups
> (`mtil`, `if(parity)` guard, inlined `std::to_string`) were applied. This file is closed.

## Goal

Make `ConservedCurrent` a first-class `LinOp` so the conserved-current kernel $K$ can be
wrapped in a `MatPoly` and applied through `MatPoly::from_cpu`, exactly like the Dirac
operators (`op_DH`, `op_Dm`, ...). This lets the production `.cu`
(`jj_conn_tpproj_claude.cu`) drop every raw `CuC*` device pointer (`d_eta`, `d_phi0`,
`d_acc`, `d_phimm`, `d_va`) and work entirely in `FermionVector` + `MatPoly` idioms.

Readability is the only objective; the numerical result must be identical.

## Design

`LinOp` (includes/sparse_matrix.h:4) interface:
```
virtual void operator()(CuC* d_res, const CuC* d_v) const = 0;
virtual void Async(CuC* d_res, const CuC* d_v, const cudaStream_t stream) const = 0;
```
`operator()` is fixed `(d_res, d_v)` -- it carries no gauge / no link. So
`ConservedCurrent` must store the gauge $U$, the link element, and the dagger flag as
state, set just before the operator is applied.

### Header changes (includes/conserved_current_claude.h)

1. `template<typename OverlapOp>` -> `template<typename OverlapOp, typename Gauge>` and
   `struct ConservedCurrent : public LinOp`. (The kernel needs $U$, whose concrete type is
   `GaugeExt<...>`; templating on it is the minimal way to store `const Gauge&`.)
2. Drop the per-method `typename Gauge` template parameter on `apply_k` / `apply_k_dag` /
   `build_W` / `build_W_wz` / the `_wz` wrappers; use the class `Gauge` type directly
   (avoids shadowing the new class template parameter). Keep `typename LinkEl` where
   present. Make all of these methods `const` (they write only through scratch pointer
   members, never reassign them -- valid in a const method).
3. Add configuration state + setters + the `LinOp` overrides:
   ```
   enum class Kind { TEMPORAL, SPATIAL };
   const Gauge* cfg_U = nullptr;
   Kind  cfg_kind = Kind::TEMPORAL;
   bool  cfg_dag  = false;
   int   cfg_ts   = 0;        // timeslice t (temporal) or spatial slot s (spatial)
   Idx   cfg_n    = 0;        // dual site n (temporal)
   BaseLink cfg_lk = {0,0};   // spatial link (spatial)

   void set_temporal(const Gauge& U, int t, Idx n, bool dag);
   void set_spatial (const Gauge& U, int s, BaseLink lk, bool dag);

   void operator()(CuC* d_res, const CuC* d_v) const override; // dispatch on cfg_kind/cfg_dag
   void Async(CuC* d_res, const CuC* d_v, const cudaStream_t) const override; // = operator() + sync
   ```
   `operator()` rebuilds the `std::pair` of the right type from the stored fields and calls
   `apply_k` / `apply_k_dag`. `Async` ignores the stream (apply_k is internally synchronous
   on the default stream) -- same posture as `LinOpDHDWrapper`.

`apply_k` / `apply_k_dag` bodies are otherwise UNCHANGED (raw-pointer device API stays
available for callers that still want it, e.g. existing micro-tests).

### Production .cu changes (jj_conn_tpproj_claude.cu)

- `ConservedCurrent<Fermion> kop(Dm);` -> `ConservedCurrent<Fermion, Gauge> kop(Dm);`
- Add `MatPoly op_K; op_K.push_back(cplx(1.0), {&kop});`.
- Delete the `CuC *d_eta,*d_phi0,*d_acc,*d_phimm,*d_va;` block, all `cudaMalloc`/`cudaFree`
  and `cudaMemcpy` for them.
- Estimator becomes (per hit), all host `FermionVector` + `from_cpu`:
  - `phi0 = D_m^{-1} eta`  : `op_DH.from_cpu(DH_eta,eta); op_DHD.solve(phi0,DH_eta,...)`
  - per site n:
    - `psi = D_m^{-dag} K^dag(n,0) eta`:
      `kop.set_temporal(U,0,n,true); op_K.from_cpu(rho_host, eta);`
      `op_Dm.from_cpu(Dm_rho, rho_host); op_DmDH.solve(psi, Dm_rho, ...)`
    - per t: `kop.set_temporal(U,t,n,false); op_K.from_cpu(phit, phi0); C[t]+=w_tp[n]*psi.dag(phit)`
  - parity (--): `phimm = tilde^{-dag} eta`; per site `v = tilde^{-dag} K^dag(n,0) phimm`;
    per t `kop.set_temporal(U,t,n,true); op_K.from_cpu(phit, va); Cmm[t]+=w_tp[n]*eta.dag(phit)`.
  - non-parity (--): `mm = conj(pp)` unchanged.

### Check-file changes (check_conserved_current_claude.cu, check_conserved_current_dag_claude.cu)

- Update instantiation: `ConservedCurrent<Fermion>` -> `ConservedCurrent<Fermion, Gauge>`
  (the `Gauge` alias already exists in both).
- Existing direct `apply_k` / `apply_k_wz` / `build_W` calls keep compiling unchanged.
- SCOPE DECISION (RESOLVED 2026-06-04): mirror the production path. The overlap-current
  Ward/conservation sections of the checks are rewritten to the `op_K.from_cpu` pattern
  (host `FermionVector`), so the tests exercise the exact production code path. Pure
  kernel-algebra micro-tests that have no production analogue (raw `build_W`, single
  `apply_k` algebra identities) may stay on the raw-pointer API.

## Chunks

1. [DONE] Header refactor: `ConservedCurrent<OverlapOp,Gauge> : public LinOp`; `set_temporal` /
   `set_spatial` / `set(U,el,dag)` + `operator()` / `Async`; kernel methods made `const`,
   per-method `Gauge` template dropped; added `<utility>` / `<cassert>`.
   Files: includes/conserved_current_claude.h
2. [DONE] Production .cu: `ConservedCurrent<Fermion,Gauge>` + `MatPoly op_K`; estimator uses
   `kop.set_temporal` + `op_K.from_cpu` + `FermionVector.dag`; all `CuC*` (d_eta/d_phi0/d_acc/
   d_phimm/d_va) and their malloc/memcpy/free removed.
   Files: jj_conn_tpproj_claude.cu
3. [DONE] Check files (mirror production path): instantiation -> `<Fermion,Gauge>`; the Ward /
   trace sections (forward tr(K) + Ward 5b/6; dag tr(K^dag) D4 + Ward D5 + adjoint D2) rewritten
   to `kop.set(...)` + `op_K.from_cpu`, dropping their input device buffers (d_ei, local d_eta,
   d_phi_dev, d_phi6, d_src, d_kout, d_psi, d_k5; D2's d_eta/d_Kphi/d_Kdeta). The pure
   kernel-algebra micro-tests (force, [X,Theta], [D_ov,Theta], D3 commutator, raw build_W) keep
   the raw-pointer API and their d_xi / d_k_xi / d_phi buffers.
   Files: check_conserved_current_claude.cu, check_conserved_current_dag_claude.cu

Sync posture (per user): the valence path uses only synchronous `from_cpu` / `on_gpu` / `solve`;
`operator()` calls the internally-synchronous `apply_k`, and `Async` adds an explicit
`cudaDeviceSynchronize`. No async call is left without an explicit sync.

User compiles/reruns the checks; Claude does not compile or run.

## Open questions

- (none open; check-file scope resolved -> mirror production path.)
