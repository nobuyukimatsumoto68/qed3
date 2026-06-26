#pragma once

#include <cmath>
#include <vector>


struct Zolotarev{
  using FLOAT = double;
  static constexpr FLOAT ONE = 1.0;
  static constexpr FLOAT TWO = 2.0;
  static constexpr FLOAT HALF = ONE/TWO;

  const int n;
  double k;
  double kp;
  const int size;

  std::vector<double> c;
  std::vector<double> cp;
  double M;
  double lambda_inv;

  double E;
  double C;
  std::vector<double> A;

  Zolotarev(const double k_=0.01,
	    const int n_=21 )
    : n(n_)
    , k(k_)
    , kp(std::sqrt(1.0-k*k))
    , size( int(n/2)+1 )
    , c(size, 0.0)
    , cp(size, 0.0)
    , A(size, 0.0)
  {
    get_coeffs();
    partial_fraction();
    C = E * 2.0 / (1.0+lambda_inv) / (k*M);
  }

  ~Zolotarev(){
  }

  void update(const double k_)
  {
    k = k_;
    kp = std::sqrt(1.0-k*k);
    get_coeffs();
    partial_fraction();
    C = E * 2.0 / (1.0+lambda_inv) / (k*M);
  }


  // A.D.Kennedy 2004; https://arxiv.org/abs/hep-lat/0402038
  static void sncndnK(const FLOAT z, const FLOAT k,
		      FLOAT& sn, FLOAT& cn, FLOAT& dn,
		      FLOAT& K) {
    FLOAT agm;
    int sgn;
    sn = arithgeom(z, ONE, std::sqrt(ONE - k*k), &agm);
    K = M_PI / (TWO * agm);
    sgn = ((int) (std::abs(z) / K)) % 4; /* sgn = 0, 1, 2, 3 */
    sgn ^= sgn >> 1; /* (sgn & 1) = 0, 1, 1, 0 */
    sgn = 1 - ((sgn & 1) << 1); /* sgn = 1, -1, -1, 1 */
    cn = ((FLOAT) sgn) * std::sqrt(ONE - sn * sn);
    dn = std::sqrt(ONE - k*k* sn * sn);
  }

  static FLOAT arithgeom(const FLOAT z, FLOAT a, FLOAT b, FLOAT* agm) {
    static FLOAT pb = -ONE;
    FLOAT xi;
    if (b <= pb) {
      pb = -ONE;
      *agm = a;
      return std::sin(z * a);
    }
    pb = b;
    xi = arithgeom(z, HALF*(a+b), std::sqrt(a*b), agm);
    return 2*a*xi / ((a+b) + (a-b)*xi*xi);
  }

  void get_coeffs() {
    double Kp = 1.0, xibar;
    for(int m=0; m<size; m++){
      double sn, cn, dn;
      double z = 2.0 * Kp * m / n;
      sncndnK( z, kp, sn, cn, dn, Kp );
      if(m==0) continue;
      c[m] = - std::pow( cn / sn, 2 );
      z = 2.0 * Kp * (m-0.5) / n;
      sncndnK( z, kp, sn, cn, dn, Kp );
      cp[m] = - std::pow( cn / sn, 2 );
      if(m==1) xibar = 1.0/dn;
    }

    M = 1.0;
    for(int m=1; m<size; m++) M *= (1.0-c[m]) / (1.0-cp[m]);

    lambda_inv = xibar / M;
    for(int m=1; m<size; m++) lambda_inv *= (1.0-c[m]*xibar*xibar) / (1.0-cp[m]*xibar*xibar);
  }

  void partial_fraction() {
    E = 1.0;
    for(int m=1; m<size; m++) {
      double numer = 1.0, denom = 1.0;
      for(int ell=1; ell<size; ell++) {
	numer *= k*k/cp[m] - k*k/c[ell];
	if(m!=ell) denom *= k*k/cp[m] - k*k/cp[ell];
      }
      A[m] = numer/denom;
      E*=c[m]/cp[m];
    }
  }

  double operator[]( const double x ) const {
    double res = 2.0 / (1.0+lambda_inv) * x / (k*M);
    for(int m=1; m<size; m++) res *= (k*k - c[m]*x*x) / (k*k - cp[m]*x*x);
    return res;
  }

  double operator()( const double x ) const {
    double res = 1.0;
    for(int m=1; m<size; m++) res += A[m] / (x*x - k*k/cp[m]);
    res *= C * x;
    return res;
  }

  inline double Delta() const {
    return (lambda_inv - 1.0) / (lambda_inv + 1.0);
  }
};


template<typename WilsonDirac>
struct OverlapWMass : public Zolotarev {
  static constexpr Idx N = Comp::N;
  static constexpr int nstreams = Comp::NSTREAMS;

  const WilsonDirac& DW;
  DWDevice<WilsonDirac> d_DW; // actual data used in M_DW, M_DWH
  CSR<N> M_DW;
  CSR<N> M_DWH;
  double lambda_max, lambda_min;

  bool is_update;
  Complex mass;

  // _claude: measure-weighted diagonal mass (mass_measure_factor_impl_plan_claude.md).
  // `mass` is now the PHYSICAL m (R=1 units; COMPLEX allowed). m_L = mass_coeff * M_mass:
  //   M_mass     = volume_matrix(pow=1) = diag(A_y / Abar)   (real, gauge-independent),
  //   mass_coeff = m * Abar / abar_s = mass * mean_dual_area / mean_ell   (complex).
  // apply_mL applies m_L; apply_mLdag applies m_L^\dagger (conj(mass_coeff)). d_mtmp = scratch.
  COO<N> M_mass;
  Complex mass_coeff = 0.0;
  CuC* d_mtmp = nullptr;

  std::vector<cudaStream_t> stream;
  std::vector<cublasHandle_t> handle;

  std::vector<CuC*> d_Ys, d_Zs, d_XYs, d_XZs;
  std::vector<CuC*> d_Ms;   // _claude: per-pole scratch for the M_mass-weighted force dot (diagonal-mass grad)
  // L1 (HMC force, hmc_force_opt_impl_plan_claude.md): X Z_m / X Y_m are LINK-INDEPENDENT, so precompute
  // them ONCE per force eval in precalc_grad (here) and reuse in grad_..._l1 instead of recomputing per
  // link.  X = M_DW/lambda_max.  Indexed 1..size-1 (m=0 unused).
  std::vector<CuC*> d_XZpre, d_XYpre;
  // _claude (diagonal-mass FORCE fix, grad_diag_mass_force_bug_claude.md Sec 5'): the exact massive
  // force carries the bra (1+m_L)eta THROUGH the resolvent. d_eta_bra = (1+m_L)eta (m=0 bra); for the
  // pole bra, d_Ys[m] is OVERWRITTEN in precalc to hat Y_m = R_m X^dag(1+m_L)eta = Y_m + mass_coeff*W_m,
  // W_m = R_m X^dag(M_mass eta) (ONE extra multishift solve per pole, massive only). d_Mmeta = M_mass eta;
  // d_Ws[istream] = per-stream W_m scratch. Massless (mass=0): untouched (d_eta_bra=eta, d_Ys[m]=Y_m).
  CuC* d_eta_bra = nullptr;
  CuC* d_Mmeta   = nullptr;
  std::vector<CuC*> d_Ws;
  // L2 (block grad over poles): contiguous N*(size-1) blocks (Z/Y inputs, X-precompute, coo outputs CY/CZ)
  // + device per-pole dot results. Filled in precalc_grad from the L1 buffers; consumed by grad_..._l2.
  CuC *d_Zg=nullptr, *d_Yg=nullptr, *d_XZg=nullptr, *d_XYg=nullptr, *d_CY=nullptr, *d_CZ=nullptr;
  CuC *d_dotA=nullptr, *d_dotB=nullptr;
  // _claude (diagonal-mass force fold for grad_l2/l4): M_mass-weighted blocks + their per-pole dots.
  // d_MA = M_mass (X Z_m) block, d_MB = M_mass (coo Z_m) block; d_dotAM/d_dotBM = the second block_dot.
  CuC *d_MA=nullptr, *d_MB=nullptr, *d_dotAM=nullptr, *d_dotBM=nullptr;
  // L4 (skip coo.do_it()): device buffers for the raw single-link COO entries (row/col/val); the
  // per-link grad COO is tiny (few entries) so MAXENT is a generous bound.
  static constexpr int MAXENT = 256;
  Idx *d_ent_i=nullptr, *d_ent_j=nullptr;  CuC *d_ent_v=nullptr;

  // contiguous N*(size-1) blocks for the multi-shift inner solve (column-major
  // [m*N + i], pole index m = 0..size-2). Used by mult/adj_deviceAsyncLaunch_ms.
  CuC *d_Zblock, *d_Yblock;

  bool is_precalc;

  explicit OverlapWMass( const WilsonDirac& DW_,
                          const Complex mass_,
                          const int n_=21,
                          const double k_=0.01,
                          const bool locate_on_gpu=true)
    : Zolotarev(k_, n_)
    , DW(DW_)
    , d_DW(DW)
    , stream(nstreams)
    , handle(nstreams)
    , d_Ys(size)
    , d_Zs(size)
    , d_XYs(size)
    , d_XZs(size)
    , d_XZpre(size)
    , d_XYpre(size)
    , d_Ms(size)
    , d_Ws(nstreams)
    , is_precalc(false)
    , is_update(true)
    , mass(mass_)
  {
    d_DW.associateCSR( M_DW, false );
    d_DW.associateCSR( M_DWH, true );

    for(int istream=0; istream<nstreams; istream++) {
      CUDA_CHECK(cudaStreamCreate( &stream[istream] ));
      CUBLAS_CHECK(cublasCreate(&handle[istream]));
      CUBLAS_CHECK(cublasSetStream(handle[istream], stream[istream]));
    }

    for(int m=0; m<size; m++) {
      CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
      CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
      CUDA_CHECK(cudaMalloc(&d_XZs[m], N*CD));
      CUDA_CHECK(cudaMalloc(&d_XYs[m], N*CD));
      CUDA_CHECK(cudaMalloc(&d_XZpre[m], N*CD));   // L1: X Z_m (link-indep precompute)
      CUDA_CHECK(cudaMalloc(&d_XYpre[m], N*CD));   // L1: X Y_m
      CUDA_CHECK(cudaMalloc(&d_Ms[m], N*CD));      // _claude: M_mass-weighted force-dot scratch
    }
    CUDA_CHECK(cudaMalloc(&d_Zblock, (size_t)N*(size-1)*CD));
    CUDA_CHECK(cudaMalloc(&d_Yblock, (size_t)N*(size-1)*CD));
    const size_t Ng = (size_t)N*(size-1)*CD;   // L2 contiguous blocks
    CUDA_CHECK(cudaMalloc(&d_Zg, Ng));  CUDA_CHECK(cudaMalloc(&d_Yg, Ng));
    CUDA_CHECK(cudaMalloc(&d_XZg, Ng)); CUDA_CHECK(cudaMalloc(&d_XYg, Ng));
    CUDA_CHECK(cudaMalloc(&d_CY, Ng));  CUDA_CHECK(cudaMalloc(&d_CZ, Ng));
    CUDA_CHECK(cudaMalloc(&d_dotA, (size_t)(size-1)*CD));
    CUDA_CHECK(cudaMalloc(&d_dotB, (size_t)(size-1)*CD));
    CUDA_CHECK(cudaMalloc(&d_MA, Ng));  CUDA_CHECK(cudaMalloc(&d_MB, Ng));   // _claude: M_mass-weighted blocks
    CUDA_CHECK(cudaMalloc(&d_dotAM, (size_t)(size-1)*CD));
    CUDA_CHECK(cudaMalloc(&d_dotBM, (size_t)(size-1)*CD));
    CUDA_CHECK(cudaMalloc(&d_ent_i, MAXENT*sizeof(Idx)));   // L4
    CUDA_CHECK(cudaMalloc(&d_ent_j, MAXENT*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_ent_v, MAXENT*CD));
    CUDA_CHECK(cudaMalloc(&d_eta_bra, N*CD));               // _claude: (1+m_L)eta bra-side force input
    CUDA_CHECK(cudaMalloc(&d_Mmeta,   N*CD));               // _claude: M_mass eta scratch
    for(int s=0; s<nstreams; s++) CUDA_CHECK(cudaMalloc(&d_Ws[s], N*CD));  // _claude: per-stream W_m

    // _claude: build the measure-weighted diagonal mass once (gauge-independent).
    DW.volume_matrix( M_mass.en, 1.0 );                                       // diag(A_y / Abar)
    M_mass.do_it();
    mass_coeff = mass * (DW.lattice.mean_dual_area / DW.lattice.mean_ell); // m * Abar/abar_s (complex)
    CUDA_CHECK(cudaMalloc(&d_mtmp, N*CD));                                     // scratch for apply_mL
  }

  ~OverlapWMass()
  {
    for(int m=0; m<size; m++) {
      CUDA_CHECK(cudaFree(d_Zs[m]));
      CUDA_CHECK(cudaFree(d_Ys[m]));
      CUDA_CHECK(cudaFree(d_XZs[m]));
      CUDA_CHECK(cudaFree(d_XYs[m]));
      CUDA_CHECK(cudaFree(d_XZpre[m]));   // L1
      CUDA_CHECK(cudaFree(d_XYpre[m]));   // L1
      CUDA_CHECK(cudaFree(d_Ms[m]));      // _claude
    }
    CUDA_CHECK(cudaFree(d_Zblock));
    CUDA_CHECK(cudaFree(d_Yblock));
    CUDA_CHECK(cudaFree(d_Zg)); CUDA_CHECK(cudaFree(d_Yg)); CUDA_CHECK(cudaFree(d_XZg));   // L2
    CUDA_CHECK(cudaFree(d_XYg)); CUDA_CHECK(cudaFree(d_CY)); CUDA_CHECK(cudaFree(d_CZ));
    CUDA_CHECK(cudaFree(d_dotA)); CUDA_CHECK(cudaFree(d_dotB));
    CUDA_CHECK(cudaFree(d_MA)); CUDA_CHECK(cudaFree(d_MB));   // _claude
    CUDA_CHECK(cudaFree(d_dotAM)); CUDA_CHECK(cudaFree(d_dotBM));
    CUDA_CHECK(cudaFree(d_ent_i)); CUDA_CHECK(cudaFree(d_ent_j)); CUDA_CHECK(cudaFree(d_ent_v));   // L4
    if(d_eta_bra) CUDA_CHECK(cudaFree(d_eta_bra));   // _claude
    if(d_Mmeta)   CUDA_CHECK(cudaFree(d_Mmeta));
    for(int s=0; s<nstreams; s++) if(d_Ws[s]) CUDA_CHECK(cudaFree(d_Ws[s]));
    if(d_mtmp) CUDA_CHECK(cudaFree(d_mtmp));   // _claude
    for(int istream=0; istream<nstreams; istream++) {
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      CUDA_CHECK(cudaStreamDestroy(stream[istream]));
      CUBLAS_CHECK(cublasDestroy(handle[istream]));
    }
  }

  template<typename Gauge>
  void update( const Gauge& U ) {
    d_DW.update( U );
    compute_lambda_max();
    // 2026-06-26 (per NM): adaptive Zolotarev re-fit REMOVED -- it was confusing and a latent
    // reversibility hole (update() is called per MD step, so the old `if(is_update) Zolotarev::
    // update(...)` could re-fit the window k MID-TRAJECTORY). The window k is now FIXED at
    // construction (L4: k=0.001, n=21; L2: k=0.01, n=17) and never changes during a run. Keep only
    // a PASSIVE warning when an eigenvalue falls below the fixed window edge k*lambda_max (the
    // spike predictor) -- do NOT mutate k. `is_update` is now vestigial.
#ifdef InfoDelta
    if( lambda_min/lambda_max < this->k )
      std::clog << "# WARNING: eval below Zolotarev window: ratio = "
                << lambda_min/lambda_max << "  k = " << this->k << std::endl;
#endif
    // -- old adaptive block (commented per convention) --
    // double safe = 0.1;
    // if( lambda_min/lambda_max < 0.1*this->k){
    //   if(is_update) Zolotarev::update(0.1*lambda_min/lambda_max);
    //   std::clog << "# Smaller Delta Detected : " << Delta() << std::endl;
    // }
    is_precalc = false;
  }

  void compute_lambda_max( const double TOL=1.0e-8, const int MAXITER=500 ) {
    Complex* q;
    CUDA_CHECK( cudaMallocHost( &q, N*CD ) );
    memset(q, 0, N*CD);

    q[0] = std::sqrt(1.0/2.0);
    q[1] = std::sqrt(1.0/2.0);

    MatPoly Op;
    Op.push_back ( cplx(1.0), {&M_DW, &M_DWH} );

    CuC *d_x, *d_q;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CUDA_CHECK(cudaMalloc(&d_q, N*CD));
    CUDA_CHECK(cudaMemcpy(d_q, q, N*CD, H2D));

    Complex dot;
    double norm=1.0, mu_0=1.0, mu_m1=1.0, mu_m2=1.0;
    double lambda=100.0, lambda_old=1000.0;

    for(int i=0; i<MAXITER; i++){
      Op.on_gpu<N>( d_x, d_q );
      //
      CUDA_CHECK(cudaMemcpy(d_q, d_x, N*CD, D2D)); // stream 2
      Op.dot2self<N>(&norm, d_x); // stream 1
      //
      Op.Zdscal<N>( 1.0/std::sqrt(norm), d_q );
      Op.dot<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);

      mu_m2=mu_m1;
      mu_m1=mu_0;
      mu_0=dot.real();

      const double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
      const double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
      lambda_old = lambda;
      lambda = mu_0 - a*std::pow(r,i);

      if(std::abs(lambda_old-lambda)/std::abs(lambda)<TOL) {
#ifdef IsVerbose
	std::clog << "# lambda_max estimate escaped in i = " << i << std::endl;
#endif
	break;
      }
    }

    CUDA_CHECK(cudaMemcpy(d_q, q, N*CD, H2D));
    double lambda2=100.0, lambda2_old=1000.0;

    for(int i=0; i<MAXITER; i++){
      Op.solve<N>( d_x, d_q, Comp::TOL_OUTER );
      //
      CUDA_CHECK(cudaMemcpy(d_q, d_x, N*CD, D2D)); // stream 2
      Op.dot2self<N>(&norm, d_x);
      //
      Op.Zdscal<N>( 1.0/std::sqrt(norm), d_q );

      Op.dot<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);

      mu_m2=mu_m1;
      mu_m1=mu_0;
      mu_0=dot.real();

      const double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
      const double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
      lambda2_old = lambda2;
      lambda2 = mu_0 - a*std::pow(r,i);

      if(std::abs(lambda2_old-lambda2)/std::abs(lambda2)<TOL) {
#ifdef IsVerbose
	std::clog << "# lambda_min estimate escaped in i = " << i << std::endl;
        std::clog << lambda << " " << lambda2 << std::endl;
#endif
	break;
      }
    }

    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_q));
    CUDA_CHECK(cudaFreeHost(q));

    // lambda_min = std::min( lambda_min, std::sqrt( std::abs((1.0-100*TOL)/lambda2) ));
    // lambda_max = std::max( lambda_max, std::sqrt( std::abs((1.0+100*TOL)*lambda) ));
    lambda_min = std::sqrt( std::abs((1.0-100*TOL)/lambda2) );
    lambda_max = std::sqrt( std::abs((1.0+100*TOL)*lambda) );
  }



  // _claude: d_out = m_L . d_in, m_L = mass_coeff * (A_y/Abar) (complex diagonal).
  // M_mass is diagonal, so this is a cheap pointwise multiply (negligible vs an overlap apply).
  void apply_mL(CuC* d_out, const CuC* d_in) const {
    MatPoly Op_m;
    Op_m.push_back( cplx(mass_coeff), {&M_mass} );
    Op_m.on_gpu<N>( d_out, d_in );
  }
  // _claude: d_out = m_L^\dagger . d_in = conj(mass_coeff) * (A_y/Abar) . d_in.
  void apply_mLdag(CuC* d_out, const CuC* d_in) const {
    MatPoly Op_m;
    Op_m.push_back( cplx(std::conj(mass_coeff)), {&M_mass} );
    Op_m.on_gpu<N>( d_out, d_in );
  }

  void mult_deviceAsyncLaunch(CuC* d_res, const CuC* d_xi) const {
    // can parallelize
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) // schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = omp_get_thread_num(); // m%nstreams;
      MatPoly Op(handle[istream], stream[istream], istream);
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solveAsync<N>( d_Zs[m], d_xi, Comp::TOL_INNER );
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
    }

    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_xi, N*CD, D2D)); // E(1+Z)
    // reduction
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);

    MatPoly Op;
    Op.push_back( cplx(1.0/(lambda_max)), {&M_DW} );
    Op.on_gpu<N>( d_res, d_Zs[0] );
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // D_ov
    CUDA_CHECK(cudaDeviceSynchronize());
    // _claude: + m_L v (diagonal measure-weighted mass), was scalar + m v:
    // Taxpy_gen<CuC,CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_res, cplx(mass), d_xi, d_res); // +m v
    apply_mL(d_mtmp, d_xi);                                                                // d_mtmp = m_L v
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_mtmp, d_res);      // d_res += m_L v
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  void adj_deviceAsyncLaunch(CuC* d_res, const CuC* d_xi) const {
    CUDA_CHECK(cudaMemcpy(d_Ys[0], d_xi, N*CD, D2D)); // E(1+Y)
    {
      MatPoly OpGlob;
      OpGlob.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      OpGlob.on_gpu<N>( d_Ys[0], d_xi );
    }

    CUDA_CHECK(cudaMemcpy(d_res, d_Ys[0], N*CD, D2D));

    // can parallelize
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) // schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = omp_get_thread_num(); // m%nstreams;
      MatPoly Op(handle[istream], stream[istream], istream);
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]);
      Op.push_back ( a, {} );
      Op.solveAsync<N>( d_Ys[m], d_Ys[0], Comp::TOL_INNER );
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
    }

    // reduction
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_Ys[m], d_res);
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // D^dag_ov
    CUDA_CHECK(cudaDeviceSynchronize());
    // _claude: + m_L^dag v (complex diagonal), was scalar + m^* v:
    // Taxpy_gen<CuC,CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_res, cplx(std::conj(mass)), d_xi, d_res); // +m^* v
    apply_mLdag(d_mtmp, d_xi);                                                             // d_mtmp = m_L^dag v
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_mtmp, d_res);      // d_res += m_L^dag v
    CUDA_CHECK(cudaDeviceSynchronize());
  }


  // _claude: mult/adj/DHD AND their _ms variants are now converted to the diagonal m_L.
  // HMC force ALSO converted (2026-06-19): the default grad and the block variants
  // grad_*_l1/_l2/_l4 fold the diagonal (1+M^*) INTO the dot product (no factored-out scalar):
  //   real((1+m^*) <a|b>)  ->  real(<a|b> + conj(mass_coeff) <a| M_mass b>)  (per-site m_L weight).
  // No NEW force term (dm_L/dU = 0, decision 6); the mass already woven into the force is now
  // diagonal. grad_l2/l4 reuse mult_coo_block(M_mass) (M_mass diagonal => per-site broadcast) +
  // a second block_dot (d_dotAM/d_dotBM). Pure inversions (valence/measurement) call mult/adj/DHD.
  // ----- multi-shift variants (C3) ----------------------------------------------
  // Numerically equivalent to mult_/adj_deviceAsyncLaunch but the per-pole solveAsync
  // loop is replaced by a single multi-shift CG pass (MatPoly::solve_multishift): the
  // matrix A = (1/lambda_max^2) D_W^dag D_W and RHS are shared by all poles, the shift
  // is sigma_m = -k^2/cp[m]. Pole m's solution is column (m-1) of d_Zblock / d_Yblock
  // and reduces with weight A[m], exactly as in the originals. Validate vs the
  // originals with test_overlap_multishift_claude.cu before switching call sites.
  void mult_deviceAsyncLaunch_ms(CuC* d_res, const CuC* d_xi) const {
    MatPoly Aseed;
    Aseed.push_back( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
    std::vector<double> sigma(size-1);
    for(int m=1; m<size; m++) sigma[m-1] = -k*k/cp[m];
    Aseed.solve_multishift<N>( d_Zblock, d_xi, sigma.data(), size-1, Comp::TOL_INNER );

    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_xi, N*CD, D2D)); // E(1+Z)
    for(int m=1; m<size; m++)
      Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zblock + (size_t)(m-1)*N, d_Zs[0]);

    MatPoly Op;
    Op.push_back( cplx(1.0/(lambda_max)), {&M_DW} );
    Op.on_gpu<N>( d_res, d_Zs[0] );
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // D_ov
    CUDA_CHECK(cudaDeviceSynchronize());
    // _claude: + m_L v (diagonal), was scalar + m v:
    // Taxpy_gen<CuC,CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_res, cplx(mass), d_xi, d_res); // +m v
    apply_mL(d_mtmp, d_xi);                                                                // d_mtmp = m_L v
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_mtmp, d_res);      // d_res += m_L v
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  void adj_deviceAsyncLaunch_ms(CuC* d_res, const CuC* d_xi) const {
    CUDA_CHECK(cudaMemcpy(d_Ys[0], d_xi, N*CD, D2D)); // E(1+Y)
    {
      MatPoly OpGlob;
      OpGlob.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      OpGlob.on_gpu<N>( d_Ys[0], d_xi );
    }
    CUDA_CHECK(cudaMemcpy(d_res, d_Ys[0], N*CD, D2D));

    MatPoly Aseed;
    Aseed.push_back( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
    std::vector<double> sigma(size-1);
    for(int m=1; m<size; m++) sigma[m-1] = -k*k/cp[m];
    Aseed.solve_multishift<N>( d_Yblock, d_Ys[0], sigma.data(), size-1, Comp::TOL_INNER );

    for(int m=1; m<size; m++)
      Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_Yblock + (size_t)(m-1)*N, d_res);
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // D^dag_ov
    CUDA_CHECK(cudaDeviceSynchronize());
    // _claude: + m_L^dag v (complex diagonal), was scalar + m^* v:
    // Taxpy_gen<CuC,CuC,N><<<NBlocks, NThreadsPerBlock>>>(d_res, cplx(std::conj(mass)), d_xi, d_res); // +m^* v
    apply_mLdag(d_mtmp, d_xi);                                                             // d_mtmp = m_L^dag v
    Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_mtmp, d_res);      // d_res += m_L^dag v
    CUDA_CHECK(cudaDeviceSynchronize());
  }


  void DHD_deviceAsyncLaunch( CuC* d_res, const CuC* d_xi ) const {
    // _claude: D_m^dag D_m v = (1+M^*) mult(v) + adj((1+M)v) - (2 Re(M) + |M|^2) v,
    //   M = m_L complex diagonal,  mult(v)=(D+M)v,  adj((1+M)v)=(D+M)^dag (1+M)v.
    // From (D+M)^dag(D+M) = (1+M^*)D + D^dag(1+M) + |M|^2 (GW: D^dag D = D+D^dag).
    // Reduces to the old scalar identity (1+m^*)(D+m)v+(1+m)(D+m)^dag v-(2Re m+|m|^2)v at M->m.
    // Still 1 mult + 1 adj. OLD scalar body preserved in overlap_wmass_obsolete_claude.h.
    CuC *d_mp, *d_mm, *d_w, *d_t;
    CUDA_CHECK(cudaMalloc(&d_mp, N*CD));
    CUDA_CHECK(cudaMalloc(&d_mm, N*CD));
    CUDA_CHECK(cudaMalloc(&d_w,  N*CD));
    CUDA_CHECK(cudaMalloc(&d_t,  N*CD));

    this->mult_deviceAsyncLaunch(d_mp, d_xi);      // d_mp = (D+M) v
    CUDA_CHECK(cudaDeviceSynchronize());

    apply_mL(d_t, d_xi);                            // d_t = m_L v   (M v; reused below)
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_w, 1.0, d_xi, d_t);    // d_w = v + M v = (1+M) v
    CUDA_CHECK(cudaDeviceSynchronize());
    this->adj_deviceAsyncLaunch(d_mm, d_w);        // d_mm = (D+M)^dag (1+M) v
    CUDA_CHECK(cudaDeviceSynchronize());

    apply_mLdag(d_w, d_mp);                         // d_w = M^* (D+M)v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, 1.0, d_mp, d_w);  // d_res = (1+M^*)(D+M)v
    CUDA_CHECK(cudaDeviceSynchronize());
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, 1.0, d_mm, d_res); // += (D+M)^dag(1+M)v
    CUDA_CHECK(cudaDeviceSynchronize());
    // subtract (2 Re M + |M|^2) v = (M + M^* + |M|^2) v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, -1.0, d_t, d_res); // -= M v
    CUDA_CHECK(cudaDeviceSynchronize());
    apply_mLdag(d_w, d_xi);                         // d_w = M^* v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, -1.0, d_w, d_res); // -= M^* v
    CUDA_CHECK(cudaDeviceSynchronize());
    apply_mLdag(d_w, d_t);                          // d_w = M^*(M v) = |M|^2 v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, -1.0, d_w, d_res); // -= |M|^2 v
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaFree(d_mp)); CUDA_CHECK(cudaFree(d_mm));
    CUDA_CHECK(cudaFree(d_w));  CUDA_CHECK(cudaFree(d_t));
  }

  void DDH_deviceAsyncLaunch( CuC* d_res, const CuC* d_xi ) const {
    this->DHD_deviceAsyncLaunch(d_res, d_xi);  // DDH = DHD for normal D
  }

  // multi-shift variants of DHD/DDH (C3b): identical to the above but route the
  // (D+m) and (D+m)^dag applies through the multi-shift mult_ms/adj_ms. Used by
  // op_Dmsq_ms in test_mrhs_claude.cu for the end-to-end A/B; production call
  // sites switch to these only after sign-off.
  void DHD_deviceAsyncLaunch_ms( CuC* d_res, const CuC* d_xi ) const {
    // _claude: same diagonal identity as DHD_deviceAsyncLaunch, routed through the _ms applies:
    //   D_m^dag D_m v = (1+M^*) mult_ms(v) + adj_ms((1+M)v) - (2 Re(M) + |M|^2) v.
    CuC *d_mp, *d_mm, *d_w, *d_t;
    CUDA_CHECK(cudaMalloc(&d_mp, N*CD));
    CUDA_CHECK(cudaMalloc(&d_mm, N*CD));
    CUDA_CHECK(cudaMalloc(&d_w,  N*CD));
    CUDA_CHECK(cudaMalloc(&d_t,  N*CD));

    this->mult_deviceAsyncLaunch_ms(d_mp, d_xi);   // d_mp = (D+M) v
    CUDA_CHECK(cudaDeviceSynchronize());

    apply_mL(d_t, d_xi);                            // d_t = M v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_w, 1.0, d_xi, d_t);    // d_w = (1+M) v
    CUDA_CHECK(cudaDeviceSynchronize());
    this->adj_deviceAsyncLaunch_ms(d_mm, d_w);     // d_mm = (D+M)^dag (1+M) v
    CUDA_CHECK(cudaDeviceSynchronize());

    apply_mLdag(d_w, d_mp);                         // d_w = M^* (D+M)v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, 1.0, d_mp, d_w);  // d_res = (1+M^*)(D+M)v
    CUDA_CHECK(cudaDeviceSynchronize());
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, 1.0, d_mm, d_res); // += (D+M)^dag(1+M)v
    CUDA_CHECK(cudaDeviceSynchronize());
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, -1.0, d_t, d_res); // -= M v
    CUDA_CHECK(cudaDeviceSynchronize());
    apply_mLdag(d_w, d_xi);                         // d_w = M^* v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, -1.0, d_w, d_res); // -= M^* v
    CUDA_CHECK(cudaDeviceSynchronize());
    apply_mLdag(d_w, d_t);                          // d_w = |M|^2 v
    Taxpy_gen<CuC,double,N><<<NBlocks,NThreadsPerBlock>>>(d_res, -1.0, d_w, d_res); // -= |M|^2 v
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaFree(d_mp)); CUDA_CHECK(cudaFree(d_mm));
    CUDA_CHECK(cudaFree(d_w));  CUDA_CHECK(cudaFree(d_t));
  }

  void DDH_deviceAsyncLaunch_ms( CuC* d_res, const CuC* d_xi ) const {
    this->DHD_deviceAsyncLaunch_ms(d_res, d_xi);  // DDH = DHD for normal D
  }

  template<typename Gauge>
  void precalc_grad_deviceAsyncLaunch( const Gauge& U, const CuC* d_eta ) {
#ifdef FORCE_MULTISHIFT
    precalc_grad_deviceAsyncLaunch_ms(U, d_eta);   // multishift force; toggle = FORCE_MULTISHIFT (def'd in the .cu)
    return;
#endif
    {
      MatPoly XH;
      XH.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      XH.on_gpu<N>(d_Ys[0], d_eta);
    }

    // can parallelize omp parallel nparallel=nstreams static assert(nparallel>=nstreams)
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) // schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = omp_get_thread_num(); // m%nstreams;
      MatPoly Op(handle[istream], stream[istream], istream);
      Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
      const CuC a = cplx(-k*k/cp[m]); Op.push_back ( a, {} );

      Op.solveAsync<N>( d_Zs[m], d_eta, Comp::TOL_INNER );
      Op.solveAsync<N>( d_Ys[m], d_Ys[0], Comp::TOL_INNER );
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
    }

    // _claude (diagonal-mass force fix, grad_diag_mass_force_bug_claude.md Sec 5'): for the MASSIVE
    // operator the exact force carries the bra (1+m_L)eta THROUGH the resolvent. Overwrite
    // d_Ys[m] -> hat Y_m = R_m X^dag(1+m_L)eta = Y_m + mass_coeff * W_m, W_m = R_m X^dag(M_mass eta)
    // (ONE extra multishift solve per pole), and set d_eta_bra = (1+m_L)eta. Massless: untouched.
    if( mass != Complex(0.0,0.0) ){
      { MatPoly Mp; Mp.push_back(cplx(1.0), {&M_mass}); Mp.on_gpu<N>(d_Mmeta, d_eta); }   // M_mass eta
      CUDA_CHECK(cudaMemcpy(d_eta_bra, d_eta, N*CD, D2D));
      Taxpy_gen<CuC,CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_eta_bra, cplx(mass_coeff), d_Mmeta, d_eta_bra); // (1+m_L)eta
      { MatPoly XH; XH.push_back(cplx(1.0/lambda_max), {&M_DWH}); XH.on_gpu<N>(d_mtmp, d_Mmeta); }   // d_mtmp = X^dag(M_mass eta)
      CUDA_CHECK(cudaDeviceSynchronize());
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams)
#endif
      for(int m=1; m<size; m++) {
        const int istream = omp_get_thread_num();
        MatPoly Op(handle[istream], stream[istream], istream);
        Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
        Op.push_back ( cplx(-k*k/cp[m]), {} );
        Op.solveAsync<N>( d_Ws[istream], d_mtmp, Comp::TOL_INNER );                    // W_m = R_m X^dag(M_mass eta)
        Taxpy_gen<CuC,CuC,N><<<NBlocks,NThreadsPerBlock,0,stream[istream]>>>(d_Ys[m], cplx(mass_coeff), d_Ws[istream], d_Ys[m]); // hat Y_m
        CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      }
    }
    else {
      CUDA_CHECK(cudaMemcpy(d_eta_bra, d_eta, N*CD, D2D));   // massless: bra = eta
    }

    // L1: precompute the LINK-INDEPENDENT X Z_m / X Y_m (X = M_DW/lambda_max) once, for grad_..._l1.
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams)
#endif
    for(int m=1; m<size; m++) {
      const int istream = omp_get_thread_num();
      MatPoly X(handle[istream], stream[istream], istream);
      X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
      X.on_gpuAsync<N>( d_XZpre[m], d_Zs[m] );
      X.on_gpuAsync<N>( d_XYpre[m], d_Ys[m] );
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
    }

    // L2: gather contiguous blocks (Z/Y inputs + X-precompute) for grad_..._l2's block-COO matvec
    for(int m=1; m<size; m++){
      CUDA_CHECK(cudaMemcpy(d_Zg  + (size_t)(m-1)*N, d_Zs[m],    N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_Yg  + (size_t)(m-1)*N, d_Ys[m],    N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_XZg + (size_t)(m-1)*N, d_XZpre[m], N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_XYg + (size_t)(m-1)*N, d_XYpre[m], N*CD, D2D));
    }

    is_precalc = true;
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  // multishift force precompute (C-force): the two SAME-RHS pole loops of precalc_grad
  // (RHS d_eta and RHS d_Ys[0]) each become ONE solve_multishift pass over the seed
  // A=(1/lambda_max^2) D_W^dag D_W with shifts -k^2/cp[m]; results copied into
  // d_Zs[1..]/d_Ys[1..] so grad_deviceAsyncLaunch is byte-identical. Reuses d_Zblock/d_Yblock.
  template<typename Gauge>
  void precalc_grad_deviceAsyncLaunch_ms( const Gauge& U, const CuC* d_eta ) {
    {
      MatPoly XH;
      XH.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
      XH.on_gpu<N>(d_Ys[0], d_eta);                  // d_Ys[0] = (1/lambda_max) D_W^dag eta
    }

    MatPoly Aseed;
    Aseed.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
    std::vector<double> sigma(size-1);
    for(int m=1; m<size; m++) sigma[m-1] = -k*k/cp[m];

    Aseed.solve_multishift<N>( d_Zblock, d_eta,   sigma.data(), size-1, Comp::TOL_INNER ); // Z_m = R_m eta
    for(int m=1; m<size; m++) CUDA_CHECK(cudaMemcpy(d_Zs[m], d_Zblock + (size_t)(m-1)*N, N*CD, D2D));

    Aseed.solve_multishift<N>( d_Yblock, d_Ys[0], sigma.data(), size-1, Comp::TOL_INNER ); // Y_m = R_m (X^dag eta)
    for(int m=1; m<size; m++) CUDA_CHECK(cudaMemcpy(d_Ys[m], d_Yblock + (size_t)(m-1)*N, N*CD, D2D));

    // _claude (diagonal-mass force fix, Sec 5'): massive -> hat Y_m = Y_m + mass_coeff R_m X^dag(M_mass eta)
    // (ONE extra multishift solve), d_eta_bra = (1+m_L)eta. Massless: d_eta_bra = eta, d_Ys[m] = Y_m.
    if( mass != Complex(0.0,0.0) ){
      { MatPoly Mp; Mp.push_back(cplx(1.0), {&M_mass}); Mp.on_gpu<N>(d_Mmeta, d_eta); }   // M_mass eta
      CUDA_CHECK(cudaMemcpy(d_eta_bra, d_eta, N*CD, D2D));
      Taxpy_gen<CuC,CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_eta_bra, cplx(mass_coeff), d_Mmeta, d_eta_bra);
      { MatPoly XH; XH.push_back(cplx(1.0/lambda_max), {&M_DWH}); XH.on_gpu<N>(d_mtmp, d_Mmeta); }   // X^dag(M_mass eta)
      Aseed.solve_multishift<N>( d_Yblock, d_mtmp, sigma.data(), size-1, Comp::TOL_INNER ); // W_m = R_m X^dag(M_mass eta)
      for(int m=1; m<size; m++)
        Taxpy_gen<CuC,CuC,N><<<NBlocks,NThreadsPerBlock>>>(d_Ys[m], cplx(mass_coeff), d_Yblock + (size_t)(m-1)*N, d_Ys[m]); // hat Y_m
      CUDA_CHECK(cudaDeviceSynchronize());
    }
    else {
      CUDA_CHECK(cudaMemcpy(d_eta_bra, d_eta, N*CD, D2D));
    }

    // L1: precompute the LINK-INDEPENDENT X Z_m / X Y_m once, for grad_..._l1 (same as the pole-loop variant).
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams)
#endif
    for(int m=1; m<size; m++) {
      const int istream = omp_get_thread_num();
      MatPoly X(handle[istream], stream[istream], istream);
      X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
      X.on_gpuAsync<N>( d_XZpre[m], d_Zs[m] );
      X.on_gpuAsync<N>( d_XYpre[m], d_Ys[m] );
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
    }

    // L2: gather contiguous blocks (Z/Y inputs + X-precompute) for grad_..._l2's block-COO matvec
    for(int m=1; m<size; m++){
      CUDA_CHECK(cudaMemcpy(d_Zg  + (size_t)(m-1)*N, d_Zs[m],    N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_Yg  + (size_t)(m-1)*N, d_Ys[m],    N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_XZg + (size_t)(m-1)*N, d_XZpre[m], N*CD, D2D));
      CUDA_CHECK(cudaMemcpy(d_XYg + (size_t)(m-1)*N, d_XYpre[m], N*CD, D2D));
    }

    is_precalc = true;
    CUDA_CHECK(cudaDeviceSynchronize());
  }

  template<typename Link, typename Gauge>
  double grad_deviceAsyncLaunch( const Link& link, const Gauge& U, const CuC* d_eta ) const {
#ifdef GRAD_L4
    return grad_deviceAsyncLaunch_l4(link, U, d_eta);   // L4: L2 + skip do_it (single-link matvec); -DGRAD_L4
#endif
#ifdef GRAD_L2
    return grad_deviceAsyncLaunch_l2(link, U, d_eta);   // L2: block over poles (precedence over L1); -DGRAD_L2
#endif
#ifdef GRAD_L1
    return grad_deviceAsyncLaunch_l1(link, U, d_eta);   // L1: hoisted X Z_m/X Y_m; toggle = -DGRAD_L1
#endif
    assert( is_precalc );

    COO<N> coo;
    DW.d_coo_format(coo.en, U, link);
    coo.do_it();

    // _claude (diagonal-mass force fix, Sec 5'): NO mass fold here. d_Ys[m] is hat Y_m =
    // R_m X^dag(1+m_L)eta (the bra carries (1+m_L) through the resolvent; = Y_m for massless), set in
    // precalc; the m=0 bra is d_eta_bra = (1+m_L)eta. Pure massless contraction of hatY (bra) vs Z (ket).
    std::vector<double> tmp2reduce(size, 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) // schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = omp_get_thread_num(); // m%nstreams;
      MatPoly X(handle[istream], stream[istream], istream);
      X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );

      X.on_gpuAsync<N>(d_XZs[m], d_Zs[m]);
      coo.Async( d_XYs[m], d_Ys[m], stream[istream] );
      CuC inner;
      X.dotAsync<N>( &inner,  d_XYs[m], d_XZs[m] );        // <coo hatY_m | X Z_m>
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      tmp2reduce[m] -= cuCreal(inner);

      X.on_gpuAsync<N>(d_XYs[m], d_Ys[m]);
      coo.Async( d_XZs[m], d_Zs[m], stream[istream] );
      X.dotAsync<N>( &inner,  d_XYs[m], d_XZs[m] );        // <X hatY_m | coo Z_m>
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      tmp2reduce[m] -= cuCreal(inner);

      tmp2reduce[m] *= A[m];
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    double res = 0.0;
    // reductions
    for(int m=1; m<size; m++) res+=tmp2reduce[m];

    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_eta, N*CD, D2D));
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);
    coo( d_Ys[0], d_Zs[0] );
    CuC inner;
    { MatPoly XH; XH.dot<N>( &inner,  d_eta_bra, d_Ys[0] ); }   // bra = (1+m_L)eta
    res += cuCreal(inner);
    res *= -2.0*C/lambda_max;
    return res;
  }

  // L1 (HMC force opt, hmc_force_opt_impl_plan_claude.md): byte-identical variant of
  // grad_deviceAsyncLaunch that DROPS the two link-independent X applies (X Z_m, X Y_m) and reads the
  // precomputed d_XZpre[m]/d_XYpre[m] (set in precalc_grad). Only the per-link COO matvecs (coo Z_m,
  // coo Y_m) remain in the pole loop. Requires precalc_grad to have run (it fills d_XZpre/d_XYpre).
  // Toggle at the call site (e.g. #define GRAD_L1 in the .cu); original grad kept as the dH reference.
  template<typename Link, typename Gauge>
  double grad_deviceAsyncLaunch_l1( const Link& link, const Gauge& U, const CuC* d_eta ) const {
    assert( is_precalc );

    COO<N> coo;
    DW.d_coo_format(coo.en, U, link);
    coo.do_it();

    // _claude (Sec 5'): no mass fold. d_Ys[m]=hatY_m, d_XYpre[m]=X hatY_m (built from hatY_m in precalc),
    // d_eta_bra=(1+m_L)eta. coo Y_m = coo hatY_m (bra); X Z_m / coo Z_m are ket (plain Z).
    std::vector<double> tmp2reduce(size, 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nstreams) // schedule(static)
#endif
    for(int m=1; m<size; m++) {
      const int istream = omp_get_thread_num(); // m%nstreams;
      MatPoly X(handle[istream], stream[istream], istream);   // for dotAsync (handle/stream only; no apply)

      coo.Async( d_XYs[m], d_Ys[m], stream[istream] );         // coo hatY_m  (link-dep)
      CuC inner;
      X.dotAsync<N>( &inner, d_XYs[m], d_XZpre[m] );           // <coo hatY_m, X Z_m>
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      tmp2reduce[m] -= cuCreal(inner);

      coo.Async( d_XZs[m], d_Zs[m], stream[istream] );         // coo Z_m  (link-dep)
      X.dotAsync<N>( &inner, d_XYpre[m], d_XZs[m] );           // <X hatY_m, coo Z_m>
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
      tmp2reduce[m] -= cuCreal(inner);

      tmp2reduce[m] *= A[m];
      CUDA_CHECK(cudaStreamSynchronize(stream[istream]));
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    double res = 0.0;
    for(int m=1; m<size; m++) res+=tmp2reduce[m];

    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_eta, N*CD, D2D));
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);
    coo( d_Ys[0], d_Zs[0] );
    CuC inner;
    { MatPoly XH; XH.dot<N>( &inner, d_eta_bra, d_Ys[0] ); }   // bra = (1+m_L)eta
    res += cuCreal(inner);
    res *= -2.0*C/lambda_max;
    return res;
  }

  // L2 (HMC force opt): BLOCK the pole loop. The per-link work after L1 is npole COO matvecs of the SAME
  // single-link CSR (coo Y_m, coo Z_m) + npole dots -- batched into 2 block-COO matvecs (mult_coo_block)
  // + 2 block-dots (block_dot), removing the npole streamed kernels + their per-pole host syncs (the
  // overhead the profile showed dominates grad). Uses the contiguous d_Zg/d_Yg/d_XZg/d_XYg from
  // precalc_grad. Numerically == grad to ~reduction-order (block_dot vs cublasZdotc); validate <1e-8.
  template<typename Link, typename Gauge>
  double grad_deviceAsyncLaunch_l2( const Link& link, const Gauge& U, const CuC* d_eta ) const {
    assert( is_precalc );
    const int npole = size-1;

    COO<N> coo;
    DW.d_coo_format(coo.en, U, link);
    coo.do_it();

    const int gblk = (int)(( (size_t)N*npole + NThreadsPerBlock - 1 ) / NThreadsPerBlock);
    // CY = coo applied to the Y-block (all poles) ; CZ = coo applied to the Z-block
    mult_coo_block<CuC,N><<<gblk,NThreadsPerBlock>>>(d_CY, d_Yg, coo.d_val, coo.d_cols, coo.d_rows, npole);
    mult_coo_block<CuC,N><<<gblk,NThreadsPerBlock>>>(d_CZ, d_Zg, coo.d_val, coo.d_cols, coo.d_rows, npole);
    // a_m = <coo Y_m, X Z_m> = conj(CY_m) . XZg_m ;  b_m = <X Y_m, coo Z_m> = conj(XYg_m) . CZ_m
    // a_m = <coo hatY_m, X Z_m> ; b_m = <X hatY_m, coo Z_m>. d_Yg/d_XYg are hatY (bra carries (1+m_L)
    // via precalc, Sec 5') -> NO mass fold.
    block_dot<N><<<npole,NThreadsPerBlock>>>(d_dotA, d_CY, d_XZg, npole);
    block_dot<N><<<npole,NThreadsPerBlock>>>(d_dotB, d_XYg, d_CZ, npole);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<CuC> hA(npole), hB(npole);
    CUDA_CHECK(cudaMemcpy(hA.data(),  d_dotA,  (size_t)npole*CD, D2H));
    CUDA_CHECK(cudaMemcpy(hB.data(),  d_dotB,  (size_t)npole*CD, D2H));

    double res = 0.0;
    for(int m=1; m<size; m++){
      res += A[m]*( -cuCreal(hA[m-1]) - cuCreal(hB[m-1]) );
    }

    // post-loop (m=0 term): bra = d_eta_bra = (1+m_L)eta
    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_eta, N*CD, D2D));
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);
    coo( d_Ys[0], d_Zs[0] );
    CuC inner;
    { MatPoly XH; XH.dot<N>( &inner, d_eta_bra, d_Ys[0] ); }
    res += cuCreal(inner);
    res *= -2.0*C/lambda_max;
    return res;
  }

  // L4 (HMC force opt): grad_l2 + SKIP coo.do_it(). Upload the few raw single-link COO entries and apply
  // them with link_matvec_block (no CSR build => no O(N) rows loop, no 3 cudaMalloc, ~11% of grad_l2).
  // Numerically == grad to ~reduction/atomic order (~1e-15); validate <1e-8.
  template<typename Link, typename Gauge>
  double grad_deviceAsyncLaunch_l4( const Link& link, const Gauge& U, const CuC* d_eta ) const {
    assert( is_precalc );
    const int npole = size-1;

    COO<N> coo;
    DW.d_coo_format(coo.en, U, link);             // host entries only; NO do_it
    const int nent = (int)coo.en.size();
    assert(nent <= MAXENT);
    std::vector<Idx> hi(nent), hj(nent);  std::vector<CuC> hv(nent);
    for(int k=0;k<nent;k++){ hi[k]=coo.en[k].i; hj[k]=coo.en[k].j; hv[k]=coo.en[k].v; }
    CUDA_CHECK(cudaMemcpy(d_ent_i, hi.data(), (size_t)nent*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_ent_j, hj.data(), (size_t)nent*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_ent_v, hv.data(), (size_t)nent*CD,          H2D));

    const int gblk = (int)(( (size_t)nent*npole + NThreadsPerBlock - 1 )/NThreadsPerBlock);
    // CY = link applied to Y-block ; CZ = link applied to Z-block  (custom single-link matvec)
    CUDA_CHECK(cudaMemset(d_CY, 0, (size_t)N*npole*CD));
    CUDA_CHECK(cudaMemset(d_CZ, 0, (size_t)N*npole*CD));
    link_matvec_block<N><<<gblk,NThreadsPerBlock>>>(d_CY, d_Yg, d_ent_i, d_ent_j, d_ent_v, nent, npole);
    link_matvec_block<N><<<gblk,NThreadsPerBlock>>>(d_CZ, d_Zg, d_ent_i, d_ent_j, d_ent_v, nent, npole);
    // d_Yg/d_XYg are hatY (bra carries (1+m_L) via precalc, Sec 5') -> NO mass fold.
    block_dot<N><<<npole,NThreadsPerBlock>>>(d_dotA, d_CY, d_XZg, npole);
    block_dot<N><<<npole,NThreadsPerBlock>>>(d_dotB, d_XYg, d_CZ, npole);
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<CuC> hA(npole), hB(npole);
    CUDA_CHECK(cudaMemcpy(hA.data(),  d_dotA,  (size_t)npole*CD, D2H));
    CUDA_CHECK(cudaMemcpy(hB.data(),  d_dotB,  (size_t)npole*CD, D2H));

    double res = 0.0;
    for(int m=1; m<size; m++){
      res += A[m]*( -cuCreal(hA[m-1]) - cuCreal(hB[m-1]) );
    }

    // post-loop (m=0 term): single-col link matvec (replaces coo(d_Ys[0], d_Zs[0])); bra = (1+m_L)eta
    CUDA_CHECK(cudaMemcpy(d_Zs[0], d_eta, N*CD, D2D));
    for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);
    CUDA_CHECK(cudaMemset(d_Ys[0], 0, N*CD));
    const int gent1 = (int)(( (size_t)nent + NThreadsPerBlock - 1 )/NThreadsPerBlock);
    link_matvec_block<N><<<gent1,NThreadsPerBlock>>>(d_Ys[0], d_Zs[0], d_ent_i, d_ent_j, d_ent_v, nent, 1);
    CUDA_CHECK(cudaDeviceSynchronize());
    CuC inner;
    { MatPoly XH; XH.dot<N>( &inner, d_eta_bra, d_Ys[0] ); }
    res += cuCreal(inner);
    res *= -2.0*C/lambda_max;
    return res;
  }

};




//   void compute_lambda_max( const double TOL=1.0e-8, const int MAXITER=500 ) {
//     std::vector<Complex> q(N, 0.0);

//     q[0] = std::sqrt(1.0/2.0);
//     q[1] = std::sqrt(1.0/2.0);

//     MatPoly Op( handle[0], stream[0] );
//     Op.push_back ( cplx(1.0), {&M_DW, &M_DWH} );

//     CuC *d_x, *d_q;
//     // CUDA_CHECK(cudaMallocAsync(&d_x, N*CD, stream[0]));
//     // CUDA_CHECK(cudaMallocAsync(&d_q, N*CD, stream[nstreams-1]));
//     // CUDA_CHECK(cudaDeviceSynchronize());
//     CUDA_CHECK(cudaMalloc(&d_x, N*CD));
//     CUDA_CHECK(cudaMalloc(&d_q, N*CD));

//     // CUDA_CHECK(cudaMemcpyAsync(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D, stream[nstreams-1]));
//     CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));
//     // CUDA_CHECK(cudaDeviceSynchronize());

//     Complex dot;
//     double norm=1.0, mu_0=1.0, mu_m1=1.0, mu_m2=1.0;
//     double lambda=100.0, lambda_old=1000.0;

//     for(int i=0; i<MAXITER; i++){
//       // Op.on_gpuAsync<N>( d_x, d_q );
//       // CUDA_CHECK(cudaDeviceSynchronize());
//       Op.on_gpu<N>( d_x, d_q );
//       //
//       // CUDA_CHECK(cudaMemcpyAsync(d_q, d_x, N*CD, D2D, stream[nstreams-1])); // stream 2
//       // Op.dot2selfAsync<N>(&norm, d_x); // stream 1
//       // CUDA_CHECK(cudaDeviceSynchronize());
//       CUDA_CHECK(cudaMemcpyAsync(d_q, d_x, N*CD, D2D, stream[nstreams-1])); // stream 2
//       Op.dot2selfAsync<N>(&norm, d_x); // stream 1
//       //
//       Op.ZdscalAsync<N>( 1.0/std::sqrt(norm), d_q );
//       Op.dotAsync<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);
//       CUDA_CHECK(cudaDeviceSynchronize());

//       mu_m2=mu_m1;
//       mu_m1=mu_0;
//       mu_0=dot.real();

//       const double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
//       const double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
//       lambda_old = lambda;
//       lambda = mu_0 - a*std::pow(r,i);

//       if(std::abs(lambda_old-lambda)/std::abs(lambda)<TOL) {
// #ifdef IsVerbose
// 	std::clog << "# lambda_max estimate escaped in i = " << i << std::endl;
// #endif
// 	break;
//       }
//     }

//     CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));
//     double lambda2=100.0, lambda2_old=1000.0;

//     for(int i=0; i<MAXITER; i++){
//       Op.solveAsync<N>( d_x, d_q, Comp::TOL_OUTER );
//       CUDA_CHECK(cudaDeviceSynchronize());
//       //
//       CUDA_CHECK(cudaMemcpyAsync(d_q, d_x, N*CD, D2D, stream[nstreams-1])); // stream 2
//       Op.dot2selfAsync<N>(&norm, d_x);
//       CUDA_CHECK(cudaDeviceSynchronize());
//       //
//       Op.ZdscalAsync<N>( 1.0/std::sqrt(norm), d_q );
//       CUDA_CHECK(cudaDeviceSynchronize());

//       Op.dotAsync<N>(reinterpret_cast<CuC*>(&dot), d_x, d_q);
//       CUDA_CHECK(cudaDeviceSynchronize());

//       mu_m2=mu_m1;
//       mu_m1=mu_0;
//       mu_0=dot.real();

//       const double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
//       const double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
//       lambda2_old = lambda2;
//       lambda2 = mu_0 - a*std::pow(r,i);

//       if(std::abs(lambda2_old-lambda2)/std::abs(lambda2)<TOL) {
// #ifdef IsVerbose
// 	std::clog << "# lambda_min estimate escaped in i = " << i << std::endl;
// #endif
// 	break;
//       }
//     }

//     CUDA_CHECK(cudaFreeAsync(d_x, stream[0]));
//     CUDA_CHECK(cudaFreeAsync(d_q, stream[nstreams-1]));
//     CUDA_CHECK(cudaDeviceSynchronize());

//     // lambda_min = 0.5*std::sqrt( 1.0/lambda2 );
//     // lambda_max = 2.0*std::sqrt( lambda );
//     lambda_min = std::sqrt( (1.0-100*TOL)/lambda2 );
//     lambda_max = std::sqrt( (1.0+100*TOL)*lambda );
//     // lambda_min = 0.01; // std::sqrt( (1.0-100*TOL)/lambda2 );
//     // lambda_max = 16; // std::sqrt( (1.0+100*TOL)*lambda );
//   }




  // void mult(std::vector<Complex>& x, const std::vector<Complex>& b) const {
  //   CuC *d_x, *d_r; // , *d_tmp, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_x, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_r, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_r, reinterpret_cast<const CuC*>(b.data()), N*CD, H2D));

  //   mult_device(d_x, d_r);

  //   CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(x.data()), d_x, N*CD, D2H));
  //   CUDA_CHECK(cudaFree(d_x));
  //   CUDA_CHECK(cudaFree(d_r));
  // }


  // void mult_device(CuC* d_res, const CuC* d_xi) const {
  //   std::vector<CuC*> d_Xs(size);
  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaMalloc(&d_Xs[m], N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_Xs[0], d_xi, N*CD, D2D)); // E(1+Z)

  //   // can parallelize
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Xs[m], d_xi );
  //   }
  //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Xs[0], A[m], d_Xs[m], d_Xs[0]);

  //   MatPoly Op;
  //   // Op.Zdscal<N>(C, d_Xs[0]);

  //   Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   Op.on_gpu<N>( d_res, d_Xs[0] );
  //   // Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // 1+V
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, C, d_res, d_xi); // 1+V

  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaFree(d_Xs[m]));
  // }


  // void adj_device(CuC* d_res, const CuC* d_xi) const {
  //   std::vector<CuC*> d_Ys(size);
  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_Ys[0], d_xi, N*CD, D2D)); // E(1+Y)

  //   MatPoly OpGlob;
  //   OpGlob.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //   OpGlob.on_gpu<N>( d_Ys[0], d_xi );
  //   CUDA_CHECK(cudaMemcpy(d_res, d_Ys[0], N*CD, D2D));

  //   // can parallelize
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Ys[m], d_Ys[0] );
  //   }

  //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, A[m], d_Ys[m], d_res);
  //   OpGlob.Zdscal<N>( C, d_res );
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // A[0]=1.0

  //   for(int m=0; m<size; m++) CUDA_CHECK(cudaFree(d_Ys[m]));
  // }


    // void sq_device( CuC* d_res, const CuC* d_xi) const {
  //   CuC *d_tmp1, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

  //   this->mult_device(d_tmp1, d_xi);
  //   this->adj_device(d_tmp2, d_xi);

  //   CUDA_CHECK(cudaMemcpy(d_res, d_tmp1, N*CD, D2D));
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_tmp2, d_res); // A[0]=1.0

  //   CUDA_CHECK(cudaFree(d_tmp1));
  //   CUDA_CHECK(cudaFree(d_tmp2));
  // }



    // double grad_device( const Link& link, const Gauge& U, CuC* d_eta ) const {
  //   std::vector<CuC*> d_Ys(size), d_Zs(size);
  //   for(int m=0; m<size; m++) {
  //     CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
  //     CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
  //   }
  //   {
  //     MatPoly XH;   XH.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //     XH.on_gpu<N>(d_Ys[0], d_eta);
  //   }

  //   // can parallelize omp parallel nparallel=nstreams static assert(nparallel>=nstreams)
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]); Op.push_back ( a, {} );

  //     Op.solve<N>( d_Zs[m], d_eta );
  //     Op.solve<N>( d_Ys[m], d_Ys[0] );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   CuC inner;
  //   double res = 0.0;

  //   MatPoly X; X.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   for(int m=1; m<size; m++) {
  //     X.on_gpu<N>(d_Zs[0], d_Zs[m]);
  //     coo( d_Ys[0], d_Ys[m] );

  //     X.Zdscal<N>( A[m], d_Zs[0] );
  //     X.dot<N>( &inner, d_Zs[0], d_Ys[0] );
  //     res -= real(inner);
  //   }
  //   for(int m=1; m<size; m++) {
  //     X.on_gpu<N>(d_Ys[0], d_Ys[m]);
  //     coo( d_Zs[0], d_Zs[m] );

  //     X.Zdscal<N>( A[m], d_Ys[0] );
  //     X.dot<N>( &inner, d_Ys[0], d_Zs[0] );
  //     res -= real(inner);
  //   }

  //   CUDA_CHECK(cudaMemcpy(d_Zs[0], d_eta, N*CD, D2D));
  //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_Zs[0], A[m], d_Zs[m], d_Zs[0]);
  //   coo( d_Ys[0], d_Zs[0] );
  //   X.dot<N>( &inner, d_eta, d_Ys[0] );
  //   res += real(inner);

  //   for(int m=0; m<size; m++) {
  //     CUDA_CHECK(cudaFree(d_Zs[m]));
  //     CUDA_CHECK(cudaFree(d_Ys[m]));
  //   }
  //   res *= -2.0*C/lambda_max;
  //   return res;
  // }
