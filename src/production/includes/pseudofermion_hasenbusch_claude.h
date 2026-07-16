#pragma once

// pseudofermion_hasenbusch_claude.h (Hasenbusch C2, hasenbusch_massless_impl_plan_claude.md)
//
// Mass-preconditioned pseudofermion stack for the MASSLESS overlap HMC. Reference:
//   M. Hasenbusch, Phys. Lett. B 519 (2001) 177, hep-lat/0107019 (two-mass preconditioning),
//   generalized to a multi-mass ladder (cf. Hasenbusch-Jansen, hep-lat/0211042).
//
// Splits det(D_0^dag D_0), D_0 = D_ov massless, across an auxiliary mass ladder
//   0 = m_0 < m_1 < ... < m_K,   D_i = D_ov + m_L(m_i)  (measure-weighted diagonal mass),
// telescoping into K ratio pseudofermions + 1 heaviest-frame pseudofermion (K+1 frames total):
//   det(D_0^dag D_0) = [ prod_{i=0}^{K-1} det(D_i^dag D_i)/det(D_{i+1}^dag D_{i+1}) ] det(D_K^dag D_K).
//
// MASS CONVENTION (per NM 2026-07-13, extended for massive targets):
//   * masses[0] = the PHYSICAL target mass m (frame 0 = D_ov + m_L(m)), applied via D.set_mass(m) -- so
//     m is the true physical fermion mass (m=0 -> massless target; m!=0 -> massive). This is the mass the
//     driver reports (mRe/mIm) and tags its output dir with.
//   * masses[1..K] = the AUXILIARY Hasenbusch masses = M_mass COEFFICIENTS, applied via D.set_mass_coeff:
//     frame i>0 has m_L = masses[i] * M_mass = masses[i] * diag(A_y/Abar) (algorithmic knobs on M_mass).
// set_frame_mass(i) dispatches the right setter per frame. (For m=0, set_mass(0)==set_mass_coeff(0), so
// the massless case is unchanged.) Auxiliary frames must stay heavier than frame 0 for the split to help.
//
// ONE shared OverlapWMass serves every frame via set_frame_mass(i) (update(U) is mass-independent;
// only the diagonal m_L add differs). K = dD_ov/dU is the SAME for all frames, so the ratio-frame
// cross term (Term B) is the mass-free bilinear force validated in C1.
//
// TWO-OPERATOR SPLIT-POLE (_claude, two_operator_force_impl_plan_claude.md): D = ACTION operator (n=21,
// accurate) serves heatbath + accept/reject (gen / update_eta / S). Df = FORCE operator (fewer poles n_f),
// SAME Wilson D_W + SAME frozen window, serves the MD force ONLY. The force is the gradient of a force
// action S_f,i = phi_i^dag D_{f,i+1} (D_{f,i}^dag D_{f,i})^{-1} D_{f,i+1}^dag phi_i built entirely from Df:
// chi_f,i = D_{f,i+1}^dag phi_i, eta_f,i = (D_{f,i}^dag D_{f,i})^{-1} chi_f,i (heavy: chi_f,K = phi_K).
// Exact by Metropolis correction (Duane-Kennedy-Pendleton-Roweth, Phys. Lett. B 195 (1987) 216; cruder MD
// approx: Clark-Kennedy, PRL 98 (2007) 051601, hep-lat/0608015). phi_i is drawn ONCE (action heatbath),
// SHARED. Force side (chi_f/eta_f/Op_DHD_f/set_frame_mass_force/update_eta_force*) mirrors the action side.
//
// Per frame (ket eta_i, action S_i, heatbath, force):
//   * Heaviest i=K (standard PF for D_K):
//       phi_K = D_K^dag xi_K,   eta_K = (D_K^dag D_K)^{-1} phi_K,   S_K = phi_K^dag eta_K.
//       Force = Term A only: precalc_grad(D_K, eta_K) + compute.
//   * Ratio i=0..K-1 (light D_i preconditioned by heavier D_{i+1}):
//       S_i = phi_i^dag D_{i+1} (D_i^dag D_i)^{-1} D_{i+1}^dag phi_i.
//       chi_i = D_{i+1}^dag phi_i,  eta_i = (D_i^dag D_i)^{-1} chi_i,  S_i = chi_i^dag eta_i.
//       Heatbath: phi_i = D_{i+1}^{-dag} D_i^dag xi_i => S_i = xi_i^dag xi_i at generation, with
//         D_{i+1}^{-dag} y = D_{i+1} (D_{i+1}^dag D_{i+1})^{-1} y  (DHD-solve then mult, both at m_{i+1}).
//       Force = Term A [precalc_grad(D_i, eta_i)+compute] + Term B [precalc_grad_bilinear(phi_i,eta_i)+compute].
//
// SIGN (C1 gate, test_hasenbusch_bilinear_l1_claude.log): compute()/grad returns the FORCE = -dS/dU,
// so Term A and Term B are BOTH forces and are ADDED directly (plain +) into dSf. No extra flip.
//
// Mass-switch discipline: set_mass is O(1) and NEVER changes mid-solve (each CG runs at one mass);
// the mass only switches between the discrete D_i / D_{i+1} phases of heatbath / eta-solve / force.

#include <cmath>
#include <vector>

template<typename Fermion, typename Force>
struct HasenbuschPF {
  using T = CuC;

  static constexpr Idx N = Comp::N;

  Fermion& D;                        // ACTION operator (n=21); heatbath + accept/reject; mutated via set_mass
  Fermion& Df;                       // FORCE operator (n_f poles); MD force only. Same D_W + frozen window.
  LinOpDHDWrapper<Fermion> M_DHD;
  MatPoly Op_DHD;
  LinOpDHDWrapper<Fermion> M_DHD_f;  // _claude: force-op DHD for the FORCE eta-solve (n_f poles)
  MatPoly Op_DHD_f;

  std::vector<Complex> masses;       // {[0]=PHYSICAL target mass (set_mass), [1..K]=auxiliary M_mass COEFFS (set_mass_coeff)}
  int K;                             // number of AUXILIARY masses = masses.size()-1 (frames = K+1)

  // per-frame device buffers (index i = 0..K)
  std::vector<CuC*> d_phi;           // heatbath output phi_i (drawn ONCE with the action op; SHARED)
  std::vector<CuC*> d_eta;           // ACTION ket eta_i = (D_i^dag D_i)^{-1} chi_i; read by S() (accept/reject)
  std::vector<CuC*> d_chi;           // ACTION chi_i = D_{i+1}^dag phi_i (ratio); = phi_i for the heavy frame
  std::vector<CuC*> d_eta_f;         // FORCE ket eta_f,i = (D_{f,i}^dag D_{f,i})^{-1} chi_f,i; read by get_force_frames

  CuC *d_xi, *d_tmp, *d_w;           // heatbath scratch
  CuC *d_chif;                       // _claude: FORCE chi_f scratch (D_{f,i+1}^dag phi_i, ratio frames)
  Force ftmp;                        // per-term force accumulation scratch

  std::vector<double> xi_sqnorm;     // _claude: xi_i^dag xi_i at last gen (heatbath identity check)

  HasenbuschPF()=delete;

  // _claude: single-operator (old) call sites -> force op = action op (no split; force==action==n=21).
  template<typename Base>
  explicit HasenbuschPF( Fermion& D_, const std::vector<Complex>& masses_, Base& base )
    : HasenbuschPF( D_, D_, masses_, base ) {}

  template<typename Base>
  explicit HasenbuschPF( Fermion& D_, Fermion& Df_, const std::vector<Complex>& masses_, Base& base )
    : D(D_)
    , Df(Df_)
    , M_DHD(D)
    , M_DHD_f(Df)
    , masses(masses_)
    , K( (int)masses_.size()-1 )
    , d_phi(masses_.size())
    , d_eta(masses_.size())
    , d_chi(masses_.size())
    , d_eta_f(masses_.size())
    , ftmp(base)
    , xi_sqnorm(masses_.size(), 0.0)
  {
    assert( K >= 1 );
    // frame 0 = masses[0] = the PHYSICAL target mass (0 = massless, !=0 = massive); frames 1..K = auxiliary
    // M_mass coefficients. (No masses[0]==0 assert -- massive targets are allowed now.)
    for(int i=0; i<=K; i++){
      CUDA_CHECK(cudaMalloc(&d_phi[i], N*CD));
      CUDA_CHECK(cudaMalloc(&d_eta[i], N*CD));
      CUDA_CHECK(cudaMalloc(&d_chi[i], N*CD));
      CUDA_CHECK(cudaMalloc(&d_eta_f[i], N*CD));
    }
    CUDA_CHECK(cudaMalloc(&d_xi,  N*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
    CUDA_CHECK(cudaMalloc(&d_w,   N*CD));
    CUDA_CHECK(cudaMalloc(&d_chif, N*CD));
    Op_DHD.push_back( cplx(1.0), {&M_DHD} );
    Op_DHD_f.push_back( cplx(1.0), {&M_DHD_f} );
  }

  ~HasenbuschPF(){
    for(int i=0; i<=K; i++){
      CUDA_CHECK(cudaFree(d_phi[i]));
      CUDA_CHECK(cudaFree(d_eta[i]));
      CUDA_CHECK(cudaFree(d_chi[i]));
      CUDA_CHECK(cudaFree(d_eta_f[i]));
    }
    CUDA_CHECK(cudaFree(d_xi));
    CUDA_CHECK(cudaFree(d_tmp));
    CUDA_CHECK(cudaFree(d_w));
    CUDA_CHECK(cudaFree(d_chif));
  }

  // _claude: per-frame mass dispatch -- frame 0 = PHYSICAL target mass (set_mass); frames 1..K = auxiliary
  // M_mass COEFFICIENTS (set_mass_coeff). See the MASS CONVENTION note at the top.
  void set_frame_mass( const int i ) {
    if( i==0 ) D.set_mass( masses[0] );
    else       D.set_mass_coeff( masses[i] );
  }

  // _claude: same per-frame mass dispatch, applied to the FORCE operator Df (the MD force uses D_f).
  void set_frame_mass_force( const int i ) {
    if( i==0 ) Df.set_mass( masses[0] );
    else       Df.set_mass_coeff( masses[i] );
  }

  // heatbath: draw phi_i so that S_i = xi_i^dag xi_i at generation (per-frame check optional).
  template<class Rng>
  void gen( Rng& rng ) {
    Complex* xi;
    CUDA_CHECK( cudaMallocHost( &xi, N*CD ) );

    for(int i=0; i<=K; i++){
      rng.fill_gaussian( xi );
      CUDA_CHECK(cudaMemcpy(d_xi, xi, N*CD, H2D));

      CuC t;
      Op_DHD.dot<N>( &t, d_xi, d_xi );      // xi_i^dag xi_i for the heatbath identity S_i = xi_i^dag xi_i
      xi_sqnorm[i] = cuCreal(t);

      if(i==K){
        // heaviest frame: standard PF, phi_K = D_K^dag xi
        set_frame_mass( K );
        D.adj_deviceAsyncLaunch_ms( d_phi[K], d_xi );
      }
      else {
        // ratio frame: phi_i = D_{i+1}^{-dag} D_i^dag xi
        set_frame_mass( i );
        D.adj_deviceAsyncLaunch_ms( d_tmp, d_xi );             // d_tmp = D_i^dag xi
        set_frame_mass( i+1 );                       // D_{i+1}^{-dag} = D_{i+1} (D_{i+1}^dag D_{i+1})^{-1}
        Op_DHD.solve<N>( d_w, d_tmp, Comp::TOL_OUTER );        // d_w = (D_{i+1}^dag D_{i+1})^{-1} d_tmp
        D.mult_deviceAsyncLaunch_ms( d_phi[i], d_w );          // phi_i = D_{i+1} d_w
      }
    }
    CUDA_CHECK(cudaFreeHost(xi));

    update_eta();
  }

  // eta_i = (D_i^dag D_i)^{-1} chi_i, chi_i = D_{i+1}^dag phi_i (ratio) / phi_i (heavy).  ONE frame.
  // (C6: per-frame so a multilevel integrator can re-solve only the frames of the group it is kicking.)
  void update_eta_frame( const int i ) {
    if(i==K){
      CUDA_CHECK(cudaMemcpy(d_chi[K], d_phi[K], N*CD, D2D));  // chi_K = phi_K
      set_frame_mass( K );
      Op_DHD.solve<N>( d_eta[K], d_chi[K], Comp::TOL_OUTER );
    }
    else {
      set_frame_mass( i+1 );
      D.adj_deviceAsyncLaunch_ms( d_chi[i], d_phi[i] );       // chi_i = D_{i+1}^dag phi_i
      set_frame_mass( i );
      Op_DHD.solve<N>( d_eta[i], d_chi[i], Comp::TOL_OUTER );  // eta_i = (D_i^dag D_i)^{-1} chi_i
    }
  }
  void update_eta_frames( const int i_lo, const int i_hi ) { for(int i=i_lo; i<=i_hi; i++) update_eta_frame(i); }
  void update_eta() { update_eta_frames(0, K); }

  // _claude (two-operator split): FORCE eta_f,i = (D_{f,i}^dag D_{f,i})^{-1} chi_f,i, chi_f,i =
  // D_{f,i+1}^dag phi_i (ratio) / phi_i (heavy). Byte mirror of update_eta_frame but every op = Df, writing
  // the SEPARATE (globally allocated) d_eta_f[i]. phi_i SHARED (action heatbath); chi_f uses the d_chif
  // scratch (the force never reads chi_f). The ACTION d_eta/d_chi are UNTOUCHED. ONE frame (per group).
  void update_eta_force_frame( const int i ) {
    if(i==K){
      set_frame_mass_force( K );
      Op_DHD_f.solve<N>( d_eta_f[K], d_phi[K], Comp::TOL_OUTER );  // eta_f,K = (D_{f,K}^dag D_{f,K})^{-1} phi_K
    }
    else {
      set_frame_mass_force( i+1 );
      Df.adj_deviceAsyncLaunch_ms( d_chif, d_phi[i] );            // chi_f,i = D_{f,i+1}^dag phi_i (scratch)
      set_frame_mass_force( i );
      Op_DHD_f.solve<N>( d_eta_f[i], d_chif, Comp::TOL_OUTER );
    }
  }
  void update_eta_force_frames( const int i_lo, const int i_hi ) { for(int i=i_lo; i<=i_hi; i++) update_eta_force_frame(i); }

  // S = sum_i chi_i^dag eta_i  (plain inner product; Op_DHD.dot is cublasZdotc, no operator apply).
  double S() const {
    double res = 0.0;
    for(int i=0; i<=K; i++){
      CuC tmp;
      Op_DHD.dot<N>( &tmp, d_chi[i], d_eta[i] );
      res += cuCreal(tmp);
    }
    return res;
  }

  // per-frame action at GENERATION (should equal xi_i^dag xi_i = S_i); for the heatbath check.
  double S_frame( const int i ) const {
    CuC tmp;
    Op_DHD.dot<N>( &tmp, d_chi[i], d_eta[i] );
    return cuCreal(tmp);
  }

  // Fermion force SUMMED over frames [i_lo, i_hi] (Term A on each; Term B on the ratio frames i<K).
  // (C6: the multilevel integrator calls this with a frame RANGE = one timescale group; the whole-stack
  // force is get_force = get_force_frames(0,K).) Self-contained: own set_mass + precalc per frame.
  // _claude (two-operator split): every operator apply here is the FORCE op Df with the FORCE eta_f
  // (update_eta_force_frames must have run at the current U). This is -dS_f/dU; accept/reject uses the
  // separate accurate action S (n=21). Otherwise byte-for-byte the original (D->Df, d_eta->d_eta_f).
  // (get_force() below is the self-contained whole-stack wrapper for single-shot callers.)
  void get_force_frames( Force& dSf, const Force& U, const int i_lo, const int i_hi ) {
    bool first = true;
    for(int i=i_lo; i<=i_hi; i++){
      // Term A: operator D_{f,i}, vector eta_f,i (validated massive DHD force)
      set_frame_mass_force( i );
      Df.precalc_grad_deviceAsyncLaunch( U, d_eta_f[i] );
      if(first){ dSf.compute( U, d_eta_f[i], Df ); first = false; }
      else { ftmp.compute( U, d_eta_f[i], Df ); dSf += ftmp; }

      // Term B: ratio frames only. Bilinear bra = phi_i, ket = eta_f,i; mass-free (K = dD_ov/dU).
      // SIGN (pinned by the C2 FD gate, test_hasenbusch_force_l1_claude.log): the numerator variation
      // gives dS_B/dU = +2Re<phi|dD_ov|eta> = +d[2Re<phi|D_ov|eta>], but grad_bilinear returns
      // -d[2Re<phi|D_ov|eta>] (C1). Term A (inverse variation) carries the extra minus of dM^{-1}, so
      // grad(D_i,eta_i) = +dS_A/dU. Thus Term B enters with the OPPOSITE relative sign -> SUBTRACT.
      if(i<K){
        Df.precalc_grad_bilinear_deviceAsyncLaunch_ms( U, d_phi[i], d_eta_f[i] );
        ftmp.compute( U, d_eta_f[i], Df );
        dSf += -1.0 * ftmp;                 // Term B enters with a minus (see sign note above)
      }
    }
  }

  // whole-stack force, SELF-CONTAINED: solves the FORCE eta_f (all frames) then computes -dS_f/dU. Used by
  // single-shot callers (legacy 2-level HMCHasenbusch, force-validation gate). The ML integrator instead
  // calls update_eta_force_frames + get_force_frames per timescale group (no double solve).
  void get_force( Force& dSf, const Force& U ) {
    update_eta_force_frames( 0, K );
    get_force_frames( dSf, U, 0, K );
  }
};
