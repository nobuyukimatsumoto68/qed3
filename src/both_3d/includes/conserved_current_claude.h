#pragma once

#include <array>
#include <vector>

// BaseLink: canonical oriented spatial link {site_w, site_z}.
// Matches GaugeExt::BaseLink = std::array<int,2>.
using BaseLink = std::array<int,2>;

// Exactly conserved vector current kernel K^{wz} for two-component overlap fermions.
// Reference: eq. (3.34) of qed3int_v2-2.pdf (two-component block; four-component diagonal eq. (3.33)).
// Usage: one-end trick kernel in meson_pq_wall_v2_claude.cu and disc_claude.cu,
//        replacing mult_sigma.

template<typename OverlapOp>
struct ConservedCurrent {
  static constexpr Idx N = Comp::N;

  const OverlapOp& D;          // const ref; read-only access to DW, M_DW, M_DWH,
                               //   A[], cp[], k, C, lambda_max, size, stream[], handle[], nstreams

  std::vector<CuC*> d_Zs;     // d_Zs[0]=Z_0; d_Zs[1..size-1]=Z_ell; pre-allocated scratch
  CuC* d_tmp1;
  CuC* d_tmp2;
  CuC* d_tmp3;

  explicit ConservedCurrent(const OverlapOp& D_)
    : D(D_)
    , d_Zs(D_.size)
    , d_tmp1(nullptr)
    , d_tmp2(nullptr)
    , d_tmp3(nullptr)
  {
    for(int m=0; m<D.size; m++) CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp3, N*CD));
  }

  ~ConservedCurrent() {
    for(int m=0; m<D.size; m++) CUDA_CHECK(cudaFree(d_Zs[m]));
    CUDA_CHECK(cudaFree(d_tmp1));
    CUDA_CHECK(cudaFree(d_tmp2));
    CUDA_CHECK(cudaFree(d_tmp3));
  }

  // build_W: computes W^{wz} (eq. 3.26) as a COO for link el.
  // W^{wz}_{xy} = (1/lambda_M)(delta_{xw}delta_{yz}C_{wz} - delta_{xz}delta_{yw}C_{zw})
  // where C_{xy} = -kappa_{xy} P_{xy} U_{xy}; d_coo_format = i*C at each entry.
  // multiply by +i/lambda_M: (i/lambda_M)*(i*C) = -C/lambda_M = W.
  template<typename Gauge, typename LinkEl>
  void build_W(COO<N>& coo, const Gauge& U, const LinkEl& el) const {
    D.DW.d_coo_format(coo.en, U, el);
    for(auto& e : coo.en) {
      const Complex z(cuCreal(e.v), cuCimag(e.v));
      e.v = cplx(z * Complex(0.0, 1.0/D.lambda_max));
    }
    coo.do_it();
  }

  // build_W_wz: spatial link wrapper (eq. (3.35)).
  template<typename Gauge>
  void build_W_wz(COO<N>& coo, const Gauge& U, int s, BaseLink lk) const {
    build_W(coo, U, std::pair<int,BaseLink>{s, lk});
  }

  // apply_k: compute d_result = K(el) d_xi for link element el; eq. (3.34).
  // Works for both spatial (std::pair<int,BaseLink>) and temporal (std::pair<int,Idx>) links.
  // d_tmp1/2/3 and d_Zs[] are used as scratch; not valid after return.
  template<typename Gauge, typename LinkEl>
  void apply_k(CuC* d_result, const CuC* d_xi, const Gauge& U, const LinkEl& el) {

    // --- Step 1: Z solves ---
#ifdef _OPENMP
#pragma omp parallel for num_threads(D.nstreams)
#endif
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

    // --- Step 2: build forward and adjoint COOs ---
    // do_conjugate() lacks bounds checking so the adjoint is built via do_it()
    // with entries manually conjugated and transposed.
    COO<N> coo_W, coo_WH;
    build_W(coo_W, U, el);
    D.DW.d_coo_format(coo_WH.en, U, el);
    // adjoint W^dag: W^dag_{ji} = conj(W_{ij}) = conj(i*z/lambda_M) = -i*conj(z)/lambda_M.
    // -i*conj(a+bi)/lambda_M = (-b-ia)/lambda_M, so Re=-Im/lM, Im=-Re/lM.
    for(auto& e : coo_WH.en) {
      const Complex z(cuCreal(e.v), cuCimag(e.v));
      e.v = cplx(Complex(-z.imag()/D.lambda_max, -z.real()/D.lambda_max));
      std::swap(e.i, e.j);
    }
    coo_WH.do_it();

    // --- Term A: d_result = C * W Z_0 ---
    coo_W(d_tmp1, d_Zs[0]);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemset(d_result, 0, N*CD));
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result, D.C, d_tmp1, d_result);
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- Term B: d_result -= C * sum_m A[m] * X R_m u_m ---
    // u_m = X^dag W Z_m - W^dag X Z_m
    for(int m=1; m<D.size; m++) {
      coo_W(d_tmp1, d_Zs[m]);
      CUDA_CHECK(cudaDeviceSynchronize());
      { MatPoly XH;
        XH.push_back(cplx(1.0/D.lambda_max), {&D.M_DWH});
        XH.on_gpu<N>(d_tmp2, d_tmp1); }
      { MatPoly X;
        X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
        X.on_gpu<N>(d_tmp3, d_Zs[m]); }
      coo_WH(d_tmp1, d_tmp3);
      CUDA_CHECK(cudaDeviceSynchronize());
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_tmp2, -1.0, d_tmp1, d_tmp2);
      CUDA_CHECK(cudaDeviceSynchronize());
      { MatPoly Op;
        Op.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
        Op.push_back(cplx(-D.k*D.k/D.cp[m]), {});
        Op.solve<N>(d_tmp1, d_tmp2, Comp::TOL_INNER); }
      { MatPoly X;
        X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
        X.on_gpu<N>(d_tmp3, d_tmp1); }
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result, -D.C*D.A[m], d_tmp3, d_result);
      CUDA_CHECK(cudaDeviceSynchronize());
    }
  }

  // apply_k_wz: spatial link wrapper; eq. (3.34).
  template<typename Gauge>
  void apply_k_wz(CuC* d_result, const CuC* d_xi, const Gauge& U, int s, BaseLink lk) {
    apply_k(d_result, d_xi, U, std::pair<int,BaseLink>{s, lk});
  }

  // apply_k_dag: compute d_result = K(el)^dag d_xi for link element el.
  // Adjoint of apply_k. With X = D_W/lambda_max, R_m = (X X^dag - k^2/cp[m])^{-1} (Hermitian),
  // C = D.C, A_m = D.A[m] (real), S = I + sum_m A_m R_m (Hermitian):
  //   K^dag = C S W^dag - C sum_m A_m R_m (W^dag X - X^dag W) R_m X^dag.
  // This is a reordering of the apply_k primitives; the W / W^dag COOs and the
  // R_m / X / X^dag operators are built exactly as in apply_k. The Step-1 seed is
  // Y_m = R_m X^dag xi (instead of Z_m = R_m xi). d_Zs[1..size-1] hold Y_m;
  // d_Zs[0] is reused as temporary storage for w0 = W^dag xi.
  // d_tmp1/2/3 and d_Zs[] are used as scratch; not valid after return.
  template<typename Gauge, typename LinkEl>
  void apply_k_dag(CuC* d_result, const CuC* d_xi, const Gauge& U, const LinkEl& el) {

    // --- Step 1: seed Y_m = R_m X^dag xi (link-independent) ---
    // p = X^dag xi stored in d_tmp1, read-only by the parallel solve loop.
    { MatPoly XH;
      XH.push_back(cplx(1.0/D.lambda_max), {&D.M_DWH});
      XH.on_gpu<N>(d_tmp1, d_xi); }
#ifdef _OPENMP
#pragma omp parallel for num_threads(D.nstreams)
#endif
    for(int m=1; m<D.size; m++) {
      const int istream = omp_get_thread_num();
      MatPoly Op(D.handle[istream], D.stream[istream], istream);
      Op.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
      Op.push_back(cplx(-D.k*D.k/D.cp[m]), {});
      Op.solveAsync<N>(d_Zs[m], d_tmp1, Comp::TOL_INNER);
      CUDA_CHECK(cudaStreamSynchronize(D.stream[istream]));
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- Step 2: build forward and adjoint COOs (identical to apply_k) ---
    COO<N> coo_W, coo_WH;
    build_W(coo_W, U, el);
    D.DW.d_coo_format(coo_WH.en, U, el);
    for(auto& e : coo_WH.en) {
      const Complex z(cuCreal(e.v), cuCimag(e.v));
      e.v = cplx(Complex(-z.imag()/D.lambda_max, -z.real()/D.lambda_max));
      std::swap(e.i, e.j);
    }
    coo_WH.do_it();

    // --- Term A^dag: d_result = C S W^dag xi = C (w0 + sum_m A_m R_m w0), w0 = W^dag xi ---
    coo_WH(d_tmp1, d_xi);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemset(d_result, 0, N*CD));
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result, D.C, d_tmp1, d_result);
    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_tmp1, N*CD, D2D)); // d_Zs[0] = w0 (stable rhs for solves)
    CUDA_CHECK(cudaDeviceSynchronize());
    for(int m=1; m<D.size; m++) {
      { MatPoly Op;
        Op.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
        Op.push_back(cplx(-D.k*D.k/D.cp[m]), {});
        Op.solve<N>(d_tmp2, d_Zs[0], Comp::TOL_INNER); } // tmp2 = R_m w0
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result, D.C*D.A[m], d_tmp2, d_result);
      CUDA_CHECK(cudaDeviceSynchronize());
    }

    // --- Term B^dag: d_result -= C sum_m A_m R_m (W^dag X - X^dag W) Y_m ---
    // Y_m = d_Zs[m] (untouched by Term A^dag). r_m = W^dag X Y_m - X^dag W Y_m.
    for(int m=1; m<D.size; m++) {
      { MatPoly X;
        X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
        X.on_gpu<N>(d_tmp3, d_Zs[m]); }              // tmp3 = X Y_m
      coo_WH(d_tmp1, d_tmp3);                         // tmp1 = W^dag X Y_m = a
      CUDA_CHECK(cudaDeviceSynchronize());
      coo_W(d_tmp2, d_Zs[m]);                         // tmp2 = W Y_m
      CUDA_CHECK(cudaDeviceSynchronize());
      { MatPoly XH;
        XH.push_back(cplx(1.0/D.lambda_max), {&D.M_DWH});
        XH.on_gpu<N>(d_tmp3, d_tmp2); }              // tmp3 = X^dag W Y_m = b
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_tmp1, -1.0, d_tmp3, d_tmp1); // tmp1 = a - b = r_m
      CUDA_CHECK(cudaDeviceSynchronize());
      { MatPoly Op;
        Op.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
        Op.push_back(cplx(-D.k*D.k/D.cp[m]), {});
        Op.solve<N>(d_tmp2, d_tmp1, Comp::TOL_INNER); } // tmp2 = R_m r_m = s_m
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result, -D.C*D.A[m], d_tmp2, d_result);
      CUDA_CHECK(cudaDeviceSynchronize());
    }
  }

  // apply_k_dag_wz: spatial link wrapper for the adjoint kernel.
  template<typename Gauge>
  void apply_k_dag_wz(CuC* d_result, const CuC* d_xi, const Gauge& U, int s, BaseLink lk) {
    apply_k_dag(d_result, d_xi, U, std::pair<int,BaseLink>{s, lk});
  }

};
