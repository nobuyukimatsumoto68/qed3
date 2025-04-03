#pragma once

#include <cmath>

template<typename Fermion>
struct PseudoFermion {
  using T = CuC;

  static constexpr Complex I = Complex(0.0, 1.0);
  const int NS=2;

  Fermion& D;
  LinOpDHDWrapper<Fermion> M_DHD;
  MatPoly Op_DHD;

  CuC *d_phi, *d_eta;
  static constexpr Idx N = Comp::N;

  PseudoFermion()=delete;

  explicit PseudoFermion( Fermion& D_ )
    : D(D_)
    , M_DHD(D)
  {
    CUDA_CHECK(cudaMalloc(&d_phi, N*CD));
    CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
    Op_DHD.push_back ( cplx(1.0), {&M_DHD} );
  }

  ~PseudoFermion(){
    CUDA_CHECK(cudaFree(d_phi));
    CUDA_CHECK(cudaFree(d_eta));
  }

  template<class Rng>
  void gen( Rng& rng ) {
    std::vector<Complex> xi(N, 0.0);

    rng.fill_gaussian( xi );

    CuC *d_xi;
    CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
    CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));
    // f_DH( d_phi, d_xi );
    D.adj_deviceAsyncLaunch( d_phi, d_xi );
    CUDA_CHECK(cudaFree(d_xi));

    update_eta();
  }

  inline void update_eta() { Op_DHD.solve<N>( d_eta, d_phi, Comp::TOL_OUTER ); } // outer CG

  double S() const {
    CuC tmp;
    Op_DHD.dot<N>(&tmp, d_phi, d_eta);
    return real(tmp);
  }


  template<typename Gauge>
  inline void get_force( Gauge& pi, const Gauge& u ) const { pi.compute( u, d_eta, D ); }


};


