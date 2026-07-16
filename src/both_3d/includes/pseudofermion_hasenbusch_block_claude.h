#pragma once

// pseudofermion_hasenbusch_block_claude.h -- Nf-flavor BLOCK (mrhs) form of HasenbuschPF.
//
// This is the NSTACK-block (mrhs) version of HasenbuschPF (pseudofermion_hasenbusch_claude.h). It packs
// the N_f/2 = NSTACK independent Hasenbusch pseudofermion stacks into ONE N*NSTACK device block PER FRAME
// and drives every per-frame ACTION inversion (heatbath gen, update_eta, S) and the MD FORCE (Term A +
// Term B) through the validated block operators BlockedMat / BlockedForce (blocked_mat_claude.h,
// overlap_force_Nfblock_claude.h) instead of the serial per-flavor loop. The frame logic (mass ladder,
// ratio + heavy telescoping, set_frame_mass switches, Term B minus sign) is byte-for-byte the serial
// HasenbuschPF -- only the operator applies become block applies.
//
// Reference (physics): M. Hasenbusch, Phys. Lett. B 519 (2001) 177, hep-lat/0107019 (two-mass
// preconditioning), generalized to a multi-mass ladder (cf. Hasenbusch-Jansen, hep-lat/0211042).
// Reference (block algorithm): B. Jegerlehner, "Krylov space solvers for shifted linear systems",
// hep-lat/9612014 (inner multishift) + standard block CG (outer), as wrapped by BlockedMat.
//
// TWO-OPERATOR SPLIT (mirrors the serial): D = ACTION operator (n=21, accurate) serves heatbath +
// accept/reject (gen / update_eta / S) via its own block pool + BlockedMat blk_D. Df = FORCE operator
// (fewer poles n_f, SAME Wilson D_W + SAME frozen window) serves the MD force ONLY via its own block pool
// + BlockedMat blk_Df + BlockedForce bforce. The force is built from Df with the FORCE ket eta_f (block)
// as Term A [bforce.compute] + Term B [bforce.compute_bilinear, ratio frames, SUBTRACTED]. phi_i is drawn
// ONCE (action heatbath) and SHARED by both sides. See the serial header for the frame/sign derivation.
//
// MASS CONVENTION (identical to serial): masses[0] = the PHYSICAL target mass m (frame 0 = D_ov + m_L(m),
// applied via set_mass); masses[1..K] = AUXILIARY Hasenbusch masses = M_mass COEFFICIENTS (applied via
// set_mass_coeff). set_frame_mass(i) / set_frame_mass_force(i) dispatch the right setter per frame. Since
// BlockedMat / BlockedForce read D.mass / D.mass_coeff of the wrapped op, set_frame_mass* composes with
// the following block adj / mult / solve_sq / compute (mass-aware block applies).
//
// Column fill order in gen matches the serial per-flavor RNG draw order (NSTACK Gaussian columns per
// frame, in flavor order), so a given start config + RNG reproduces the serial phi columns.
//
// Include AFTER blocked_mat_claude.h and overlap_force_Nfblock_claude.h (for BlockMemPool / BlockedMat /
// BlockedForce and the block kernels), and the overlap header (for the Fermion op type).

#include <cmath>
#include <vector>

template<typename Fermion, typename Force, int NSTACK>
struct HasenbuschPFBlock {
  using T = CuC;

  static constexpr Idx N = Comp::N;

  Fermion& D;                        // ACTION operator (n=21); heatbath + accept/reject; mutated via set_mass
  Fermion& Df;                       // FORCE operator (n_f poles); MD force only. Same D_W + frozen window.

  std::vector<Complex> masses;       // {[0]=PHYSICAL target mass (set_mass), [1..K]=auxiliary M_mass COEFFS (set_mass_coeff)}
  int K;                             // number of AUXILIARY masses = masses.size()-1 (frames = K+1)

  BlockMemPool<N,NSTACK> pool_D;     // ACTION block scratch (npole = D.size-1)
  BlockMemPool<N,NSTACK> pool_Df;    // FORCE  block scratch (npole = Df.size-1)
  BlockedMat<N,NSTACK,Fermion> blk_D;   // action block adj / mult / solve_sq (gen / update_eta / S)
  BlockedMat<N,NSTACK,Fermion> blk_Df;  // force  block adj / solve_sq (update_eta_force)
  BlockedForce<N,NSTACK,Fermion> bforce; // force block Term A (compute) + Term B (compute_bilinear), wraps Df

  // per-frame device BLOCKS (index i = 0..K), each N*NSTACK: column c = flavor c of that frame.
  std::vector<CuC*> d_phi;           // heatbath output phi_i (drawn ONCE with the action op; SHARED)
  std::vector<CuC*> d_eta;           // ACTION ket eta_i = (D_i^dag D_i)^{-1} chi_i; read by S() (accept/reject)
  std::vector<CuC*> d_chi;           // ACTION chi_i = D_{i+1}^dag phi_i (ratio); = phi_i for the heavy frame
  std::vector<CuC*> d_eta_f;         // FORCE ket eta_f,i = (D_{f,i}^dag D_{f,i})^{-1} chi_f,i; read by get_force_frames

  CuC *d_xi_blk, *d_tmp_blk, *d_w_blk;  // heatbath scratch blocks (N*NSTACK each)
  CuC *d_chif_blk;                      // FORCE chi_f block scratch (D_{f,i+1}^dag phi_i, ratio frames)
  Complex* xi_host;                     // heatbath Gaussian staging (N*NSTACK host-pinned)

  Force ftmp;                        // per-term force accumulation scratch

  std::vector<double> xi_sqnorm;     // sum_c xi_{i,c}^dag xi_{i,c} at last gen (heatbath identity check)

  HasenbuschPFBlock()=delete;

  template<typename Base>
  explicit HasenbuschPFBlock( Fermion& D_, Fermion& Df_, const std::vector<Complex>& masses_, Base& base )
    : D(D_)
    , Df(Df_)
    , masses(masses_)
    , K( (int)masses_.size()-1 )
    , pool_D(D_.size-1, true)
    , pool_Df(Df_.size-1, true)
    , blk_D(D_, pool_D)
    , blk_Df(Df_, pool_Df)
    , bforce(Df_, blk_Df)
    , d_phi(masses_.size())
    , d_eta(masses_.size())
    , d_chi(masses_.size())
    , d_eta_f(masses_.size())
    , ftmp(base)
    , xi_sqnorm(masses_.size(), 0.0)
  {
    assert( K >= 1 );
    // frame 0 = masses[0] = the PHYSICAL target mass (0 = massless, !=0 = massive); frames 1..K = auxiliary
    // M_mass coefficients. (No masses[0]==0 assert -- massive targets are allowed.)
    for(int i=0; i<=K; i++){
      CUDA_CHECK(cudaMalloc(&d_phi[i],   (size_t)N*NSTACK*CD));
      CUDA_CHECK(cudaMalloc(&d_eta[i],   (size_t)N*NSTACK*CD));
      CUDA_CHECK(cudaMalloc(&d_chi[i],   (size_t)N*NSTACK*CD));
      CUDA_CHECK(cudaMalloc(&d_eta_f[i], (size_t)N*NSTACK*CD));
    }
    CUDA_CHECK(cudaMalloc(&d_xi_blk,   (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_tmp_blk,  (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_w_blk,    (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMalloc(&d_chif_blk, (size_t)N*NSTACK*CD));
    CUDA_CHECK(cudaMallocHost(&xi_host, (size_t)N*NSTACK*CD));
  }

  ~HasenbuschPFBlock(){
    for(int i=0; i<=K; i++){
      CUDA_CHECK(cudaFree(d_phi[i]));
      CUDA_CHECK(cudaFree(d_eta[i]));
      CUDA_CHECK(cudaFree(d_chi[i]));
      CUDA_CHECK(cudaFree(d_eta_f[i]));
    }
    CUDA_CHECK(cudaFree(d_xi_blk));
    CUDA_CHECK(cudaFree(d_tmp_blk));
    CUDA_CHECK(cudaFree(d_w_blk));
    CUDA_CHECK(cudaFree(d_chif_blk));
    CUDA_CHECK(cudaFreeHost(xi_host));
  }

  HasenbuschPFBlock(const HasenbuschPFBlock&)=delete;
  HasenbuschPFBlock& operator=(const HasenbuschPFBlock&)=delete;

  // per-frame mass dispatch (identical to serial): frame 0 = PHYSICAL target mass (set_mass); frames 1..K
  // = auxiliary M_mass COEFFICIENTS (set_mass_coeff). Acts on the ACTION op D (seen by blk_D).
  void set_frame_mass( const int i ) {
    if( i==0 ) D.set_mass( masses[0] );
    else       D.set_mass_coeff( masses[i] );
  }

  // same per-frame mass dispatch on the FORCE op Df (seen by blk_Df and bforce).
  void set_frame_mass_force( const int i ) {
    if( i==0 ) Df.set_mass( masses[0] );
    else       Df.set_mass_coeff( masses[i] );
  }

  // heatbath: draw NSTACK Gaussian columns per frame so that S_i = xi_i^dag xi_i at generation. Block
  // mirror of the serial gen: heavy frame phi_K = D_K^dag xi; ratio frame phi_i = D_{i+1}^{-dag} D_i^dag xi
  // with D_{i+1}^{-dag} = D_{i+1} (D_{i+1}^dag D_{i+1})^{-1} (block solve_sq then block mult, both at m_{i+1}).
  template<class Rng>
  void gen( Rng& rng ) {
    for(int i=0; i<=K; i++){
      for(int c=0; c<NSTACK; c++) rng.fill_gaussian( xi_host + (size_t)c*N );  // flavor-order columns
      CUDA_CHECK(cudaMemcpy(d_xi_blk, reinterpret_cast<CuC*>(xi_host), (size_t)N*NSTACK*CD, H2D));

      double s = 0.0;                                   // xi_i^dag xi_i = sum_c xi_c^dag xi_c (heatbath identity)
      for(int c=0; c<NSTACK; c++){
        CuC d;
        CUBLAS_CHECK(cublasZdotc(blk_D.handle, N, d_xi_blk + (size_t)c*N, 1, d_xi_blk + (size_t)c*N, 1, &d));
        s += real(d);
      }
      xi_sqnorm[i] = s;

      if(i==K){
        // heaviest frame: standard PF, phi_K = D_K^dag xi
        set_frame_mass( K );
        blk_D.adj( d_phi[K], d_xi_blk );
      }
      else {
        // ratio frame: phi_i = D_{i+1}^{-dag} D_i^dag xi
        set_frame_mass( i );
        blk_D.adj( d_tmp_blk, d_xi_blk );                         // d_tmp = D_i^dag xi
        set_frame_mass( i+1 );                                    // D_{i+1}^{-dag} = D_{i+1} (D_{i+1}^dag D_{i+1})^{-1}
        blk_D.solve_sq( d_w_blk, d_tmp_blk, Comp::TOL_OUTER );    // d_w = (D_{i+1}^dag D_{i+1})^{-1} d_tmp
        blk_D.mult( d_phi[i], d_w_blk );                          // phi_i = D_{i+1} d_w
      }
    }

    update_eta();
  }

  // eta_i = (D_i^dag D_i)^{-1} chi_i, chi_i = D_{i+1}^dag phi_i (ratio) / phi_i (heavy). ONE frame (block).
  void update_eta_frame( const int i ) {
    if(i==K){
      CUDA_CHECK(cudaMemcpy(d_chi[K], d_phi[K], (size_t)N*NSTACK*CD, D2D));   // chi_K = phi_K
      set_frame_mass( K );
      blk_D.solve_sq( d_eta[K], d_chi[K], Comp::TOL_OUTER );
    }
    else {
      set_frame_mass( i+1 );
      blk_D.adj( d_chi[i], d_phi[i] );                            // chi_i = D_{i+1}^dag phi_i
      set_frame_mass( i );
      blk_D.solve_sq( d_eta[i], d_chi[i], Comp::TOL_OUTER );      // eta_i = (D_i^dag D_i)^{-1} chi_i
    }
  }
  void update_eta_frames( const int i_lo, const int i_hi ) { for(int i=i_lo; i<=i_hi; i++) update_eta_frame(i); }
  void update_eta() { update_eta_frames(0, K); }

  // FORCE eta_f,i = (D_{f,i}^dag D_{f,i})^{-1} chi_f,i, chi_f,i = D_{f,i+1}^dag phi_i (ratio) / phi_i
  // (heavy). Block mirror of update_eta_frame but every op = Df (blk_Df / set_frame_mass_force), writing
  // the SEPARATE d_eta_f[i]. phi_i SHARED; chi_f uses the d_chif_blk scratch. ONE frame (per group).
  void update_eta_force_frame( const int i ) {
    if(i==K){
      set_frame_mass_force( K );
      blk_Df.solve_sq( d_eta_f[K], d_phi[K], Comp::TOL_OUTER );   // eta_f,K = (D_{f,K}^dag D_{f,K})^{-1} phi_K
    }
    else {
      set_frame_mass_force( i+1 );
      blk_Df.adj( d_chif_blk, d_phi[i] );                         // chi_f,i = D_{f,i+1}^dag phi_i (scratch)
      set_frame_mass_force( i );
      blk_Df.solve_sq( d_eta_f[i], d_chif_blk, Comp::TOL_OUTER );
    }
  }
  void update_eta_force_frames( const int i_lo, const int i_hi ) { for(int i=i_lo; i<=i_hi; i++) update_eta_force_frame(i); }

  // S = sum_i sum_c Re( chi_{i,c}^dag eta_{i,c} )  (plain cublasZdotc per column, summed over frames).
  double S() const {
    double res = 0.0;
    for(int i=0; i<=K; i++){
      for(int c=0; c<NSTACK; c++){
        CuC d;
        CUBLAS_CHECK(cublasZdotc(blk_D.handle, N, d_chi[i] + (size_t)c*N, 1, d_eta[i] + (size_t)c*N, 1, &d));
        res += real(d);
      }
    }
    return res;
  }

  // per-frame action at GENERATION (should equal xi_i^dag xi_i = S_i); for the heatbath check.
  double S_frame( const int i ) const {
    double res = 0.0;
    for(int c=0; c<NSTACK; c++){
      CuC d;
      CUBLAS_CHECK(cublasZdotc(blk_D.handle, N, d_chi[i] + (size_t)c*N, 1, d_eta[i] + (size_t)c*N, 1, &d));
      res += real(d);
    }
    return res;
  }

  // Fermion force SUMMED over frames [i_lo, i_hi] (Term A on each; Term B on the ratio frames i<K). Block
  // mirror of the serial get_force_frames: every apply is the FORCE op Df (via bforce, wrapping Df) with
  // the FORCE ket eta_f (update_eta_force_frames must have run at the current U). This is -dS_f/dU;
  // accept/reject uses the separate accurate action S. Self-contained: own set_frame_mass_force per frame.
  void get_force_frames( Force& dSf, const Force& U, const int i_lo, const int i_hi ) {
    bool first = true;
    for(int i=i_lo; i<=i_hi; i++){
      // Term A: operator D_{f,i}, block ket eta_f,i (bforce.compute does precalc + all-link grad).
      set_frame_mass_force( i );
      if(first){ bforce.compute( dSf, U, d_eta_f[i] ); first = false; }
      else { bforce.compute( ftmp, U, d_eta_f[i] ); dSf += ftmp; }

      // Term B: ratio frames only. Bilinear bra = phi_i (block), ket = eta_f,i (block); mass-free (K =
      // dD_ov/dU). compute_bilinear returns the RAW positive bilinear force; Term B enters with the
      // OPPOSITE relative sign -> SUBTRACT (pinned by the serial C2 FD gate; see the serial header note).
      if(i<K){
        bforce.compute_bilinear( ftmp, U, d_phi[i], d_eta_f[i] );
        dSf += -1.0 * ftmp;
      }
    }
  }

  // whole-stack force, SELF-CONTAINED: solves the FORCE eta_f (all frames) then computes -dS_f/dU. Used by
  // single-shot callers (2-level HMC, force-validation gate). The ML integrator instead calls
  // update_eta_force_frames + get_force_frames per timescale group (no double solve).
  void get_force( Force& dSf, const Force& U ) {
    update_eta_force_frames( 0, K );
    get_force_frames( dSf, U, 0, K );
  }
};
