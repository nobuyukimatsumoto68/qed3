#pragma once

#include <vector>

// Exactly conserved vector current kernel k^{wz} for two-component overlap fermions.
// Reference: eq. 3.15 (top-left block) of qed3int-1.pdf.
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

  // build \mathcal{W}^{wz}: forward COO for spatial link (s, lk).
  // Mirrors DiracExt::d_coo_format (dirac_ext.h:368) with the I* factor removed.
  template<typename Gauge>
  void build_W_wz(COO<N>& coo, const Gauge& U, int s, BaseLink lk) const {
    D.DW.d_coo_format(coo.en, U, {s, lk});
    // d_coo_format inserts I*exp(...); remove I by multiplying each entry by -I
    for(auto& e : coo.en) {
      const Complex z(cuCreal(e.v), cuCimag(e.v));
      e.v = cplx(z * Complex(0.0, -1.0));
    }
    coo.do_it();
  }

  // apply_k_wz: compute d_result = k^{wz} d_xi for spatial link (s, lk).
  // Step 1: parallel Z solves (mirrors precalc_grad_deviceAsyncLaunch:409-418, Z only).
  // Step 2: Term A + Term B (sequential over poles; per-pole CG for Term B).
  // d_tmp1/2/3 and d_Zs[] are used as scratch; not valid after return.
  template<typename Gauge>
  void apply_k_wz(CuC* d_result, const CuC* d_xi, const Gauge& U, int s, BaseLink lk) {

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
    build_W_wz(coo_W, U, s, lk);
    D.DW.d_coo_format(coo_WH.en, U, {s, lk});
    // adjoint \mathcal{W}^{wz\dagger}: W^dag_{ji} = conj(W_{ij}) = conj(-I*z) = I*conj(z).
    // I*conj(a+bi) = b + ia, so swap Re and Im, then swap row and col indices.
    for(auto& e : coo_WH.en) {
      const Complex z(cuCreal(e.v), cuCimag(e.v));
      e.v = cplx(Complex(z.imag(), z.real()));
      std::swap(e.i, e.j);
    }
    coo_WH.do_it();

    // --- Term A: d_result = C * \mathcal{W}^{wz} Z_0 ---
    coo_W(d_tmp1, d_Zs[0]);
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemset(d_result, 0, N*CD));
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result, D.C, d_tmp1, d_result);
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- Term B: d_result -= C * sum_m A[m] * X R_m u_m ---
    // u_m = X^dag \mathcal{W}^{wz} Z_m - \mathcal{W}^{wz\dagger} X Z_m
    for(int m=1; m<D.size; m++) {
      // d_tmp1 = W Z_m
      coo_W(d_tmp1, d_Zs[m]);
      CUDA_CHECK(cudaDeviceSynchronize());
      // d_tmp2 = X^dag d_tmp1 = X^dag W Z_m
      { MatPoly XH;
        XH.push_back(cplx(1.0/D.lambda_max), {&D.M_DWH});
        XH.on_gpu<N>(d_tmp2, d_tmp1); }
      // d_tmp3 = X Z_m
      { MatPoly X;
        X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
        X.on_gpu<N>(d_tmp3, d_Zs[m]); }
      // d_tmp1 = W^dag d_tmp3 = W^dag X Z_m
      coo_WH(d_tmp1, d_tmp3);
      CUDA_CHECK(cudaDeviceSynchronize());
      // d_tmp2 = u_m = X^dag W Z_m - W^dag X Z_m
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_tmp2, -1.0, d_tmp1, d_tmp2);
      CUDA_CHECK(cudaDeviceSynchronize());
      // d_tmp1 = R_m u_m  (CG solve)
      { MatPoly Op;
        Op.push_back(cplx(1.0/(D.lambda_max*D.lambda_max)), {&D.M_DW, &D.M_DWH});
        Op.push_back(cplx(-D.k*D.k/D.cp[m]), {});
        Op.solve<N>(d_tmp1, d_tmp2, Comp::TOL_INNER); }
      // d_tmp3 = X R_m u_m
      { MatPoly X;
        X.push_back(cplx(1.0/D.lambda_max), {&D.M_DW});
        X.on_gpu<N>(d_tmp3, d_tmp1); }
      // d_result -= C * A[m] * d_tmp3
      Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_result, -D.C*D.A[m], d_tmp3, d_result);
      CUDA_CHECK(cudaDeviceSynchronize());
    }
  }

};
