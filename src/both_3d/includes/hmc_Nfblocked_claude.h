#pragma once

// hmc_Nfblocked_claude.h -- HMC2Block: a copy of HMC2 (hmc.h) that holds a PseudoFermionBlock
// (the NSTACK = Nf/2 flavors packed as one block) instead of std::vector<PseudoFermion>. H() and run()
// are identical to HMC2 except S/gen/integrate go through the block. The FORCE is still serial per
// flavor inside MinimumNorm2Block. All template params appear in the ctor args -> C++17 CTAD deduces
// them (incl. NSTACK from bpf), matching the HMC2 usage pattern.

template <typename Rng, typename Action, typename Fermion,
          typename Gauge, typename Force, typename Integrator, int NSTACK>
struct HMC2Block {
  Rng& rng;
  Action* Sg;
  Fermion* fermion;                       // on device; needs .update( Gauge& U )

  Gauge& U;
  Force& pi;
  PseudoFermionBlock<Fermion,NSTACK>& bpf;

  Integrator* integrator;

  static constexpr Idx N = Comp::N;

  HMC2Block()=delete;

  explicit HMC2Block(Rng& rng_,
                     Action* Sg_,
                     Fermion* fermion_,
                     Gauge& U_,
                     Force& pi_,
                     PseudoFermionBlock<Fermion,NSTACK>& bpf_,
                     Integrator* integrator_)
    : rng(rng_)
    , Sg(Sg_)
    , fermion(fermion_)
    , U(U_)
    , pi(pi_)
    , bpf(bpf_)
    , integrator(integrator_)
  {}

  ~HMC2Block(){}

  double H() {
    double res = pi.squared_norm();
    res *= 0.5;
    res += Sg->operator()(U);
    res += bpf.S();
    return res;
  }

  inline void integrate() {
    integrator->integrate( U, pi, Sg, fermion, bpf );
  }

  void run( double& r,
            double& dH,
            bool& is_accept,
            const bool no_reject = false ) {
    pi.gaussian( rng );

    Gauge U0( U );
    bpf.gen( rng );                        // new phi (block heat-bath) + eta at current U

    const double h0 = H();
    integrate();
    const double h1 = H();

    dH = h1-h0;
    r = std::min( 1.0, std::exp(-dH) );
    const double a = rng.uniform();
    if( a < r || no_reject ){
      is_accept=true;
    }
    else {
      is_accept=false;
      U = U0;
      fermion->update( U );                // restore the D(U) matrix to the original
    }
  }
};
