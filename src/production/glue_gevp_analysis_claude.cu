// glue_gevp_analysis_claude.cu
//
// Streaming GEVP glueball-mass analysis for the F^2 (glue_f2_claude.cu) and linear
// F_12 (glue2_msm_claude.cu) measurements. Host-only (Eigen); reads the plain-text
// F_corr.k / F.k files, accumulates per-bin jackknife sums (memory-light: it never
// holds all configs at once), and runs the SAME GEVP pipeline as glue_ylms3_claude.ipynb:
//   fold t <-> Nt-t  ->  connected (vacuum-subtracted) matrix  ->  zero_noise
//   ->  t0-metric GEVP  Delta_eff(t) = -log lambda / (dt*at).
//
// One binary handles both operators via CLI (drop_l0 + NCH):
//   F^2 :  N_FLOW=3 NCH=2 drop_l0=0  -> n_lm=16, nops=96
//   F_12:  N_FLOW=3 NCH=1 drop_l0=1  -> n_lm=15, nops=45
// operator index op = iflow*(NCH*n_lm) + ich*n_lm + ilm  (matches both .cu drivers).
//
// Usage:
//   ./a.out <dir> <N_FLOW> <NCH> <drop_l0> <at> <outfile> [nops2=8] [tcut=5] [binsize=0auto] [kmin=100]
//
// GEVP faithfully mirrors glue_ylms3_claude.ipynb (get_eigvals / zero_noise / vacuum_sub).

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <Eigen/Dense>
#include <highfive/H5File.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// OLD (per-sample pruning: the kept dimension was decided independently inside EVERY jackknife sample --
// forbidden by NM 2026-10-12; kept for reference):
// static MatrixXd inv_sqrt_sym(const MatrixXd& M, double rtol){
//   // Regularized C^{-1/2} via symmetric eigendecomposition. Directions with eigenvalue
//   // <= rtol*max (or non-positive) are PROJECTED OUT (id=0) -- the SVD pruning the
//   // glue_ylms3 notebook recommends for the near-singular multi-flow metric. Prevents
//   // the NaN from inv(sqrtm(.)) on an indefinite matrix.
//   Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
//   VectorXd d = es.eigenvalues();
//   MatrixXd V = es.eigenvectors();
//   const double dmax = d.maxCoeff();
//   VectorXd id(d.size());
//   for(int a=0; a<d.size(); a++){
//     id(a) = (d(a) > rtol*dmax && d(a) > 0.0) ? 1.0/std::sqrt(d(a)) : 0.0;
//   }
//   return V * id.asDiagonal() * V.transpose();
// }

// UNIFORM metric pruning (NM 2026-10-12): the kept dimension of every regularized C^{-1/2} is decided ONCE on
// the FULL-sample correlator (rtol rule: eigenvalues > rtol*max and > 0) and then imposed IDENTICALLY on every
// jackknife sample (the same number of largest eigen-directions).  Calls are recorded in call order.
struct KeepDims {
  bool determine = true;     // true: decide by the rtol rule and record; false: replay the recorded counts
  std::vector<int> n;        // kept dimension per inv_sqrt_sym call, in call order
  std::size_t pos = 0;
  int n_mismatch = 0;        // replay mode: calls where the per-sample rtol rule would have kept a different count
  int n_nonpos = 0;          // replay mode: kept eigenvalues <= 0 (that sample becomes NaN -> loud, never dropped)
  void rewind(){ pos = 0; }
};

static MatrixXd inv_sqrt_sym(const MatrixXd& M, double rtol, KeepDims& kd){
  Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
  VectorXd d = es.eigenvalues();      // ascending
  MatrixXd V = es.eigenvectors();
  const int nd = (int)d.size();
  const double dmax = d.maxCoeff();
  int nrtol = 0;
  for(int a=0; a<nd; a++){
    if(d(a) > rtol*dmax && d(a) > 0.0) nrtol++;
  }
  int nkeep = nrtol;
  if(kd.determine){
    kd.n.push_back(nrtol);
  }
  else{
    nkeep = kd.n.at(kd.pos);
    if(nrtol != nkeep) kd.n_mismatch++;
  }
  kd.pos++;
  VectorXd id = VectorXd::Zero(nd);
  for(int a=nd-nkeep; a<nd; a++){     // the nkeep LARGEST eigenvalues
    if(d(a) > 0.0){
      id(a) = 1.0/std::sqrt(d(a));
    }
    else{
      kd.n_nonpos++;
      id(a) = std::nan("");
    }
  }
  return V * id.asDiagonal() * V.transpose();
}

// Connected correlator preparation shared by the full sample and every jackknife sample:
// vacuum subtraction C_ij -= <O_i><O_j> (vacsub) and zero_noise of the rotationally forbidden entries.
static void prepare_conn(MatrixXd& conn, const VectorXd& m1, int nops, int Nthalf, int vacsub,
                         const std::vector<int>& opL, const std::vector<int>& opM){
  if(vacsub){
    for(int i=0; i<nops; i++){
      for(int j=0; j<nops; j++){
        const double p = m1(i)*m1(j);
        for(int t=0; t<Nthalf; t++) conn(i*nops+j, t) -= p;
      }
    }
  }
  for(int i=0; i<nops; i++){
    for(int j=0; j<nops; j++){
      const int l1=opL[i], l2=opL[j], mm1=opM[i], mm2=opM[j];
      bool z=false;
      if((l1==1 && l2==0) || (l1==0 && l2==1)) z=true;
      else if(l1==1 && l2==1 && mm1!=mm2) z=true;
      if(z){
        for(int t=0; t<Nthalf; t++) conn(i*nops+j, t) = 0.0;
      }
    }
  }
}

// GEVP on a subgroup of operators (e.g. one m-block = the n_shapes shapes of a fixed (l,m)).
// Builds the group's folded correlator matCA[tt] (ng x ng), optionally rebases at trebase, then for
// ti=1..tcut-1 solves the GEVP and writes the top `nstates_g` effmasses into out[ti-1][0..nstates_g-1]
// (out ordered ascending in state = descending lambda; last = ground = lightest).
static void gevp_group(const std::vector<int>& gops, const MatrixXd& conn, int nops, int Nthalf,
                       int trebase, int dt, int tcut, double rtol, int use_rebase, int omit_top,
                       double at, int nstates_g, std::vector<std::vector<double>>& out,
                       int t0mode, int t0fix, const std::vector<int>& offs, KeepDims& kd)
{
  const int ng = (int)gops.size();
  std::vector<MatrixXd> matCA(Nthalf);
  for(int tt=0; tt<Nthalf; tt++){
    MatrixXd Ck(ng,ng);
    for(int a=0; a<ng; a++){
      for(int b=0; b<ng; b++) Ck(a,b) = conn(gops[a]*nops + gops[b], tt);
    }
    matCA[tt] = 0.5*(Ck + Ck.transpose());
  }

  // ---- FIXED-t0 GEVP (t0mode=1): mirrors the axial rebase_vectors + rebased_effmass_fixed algebra
  //      (Blossier et al. arXiv:0902.1265; block-Hankel/GPOF Aubin-Orginos arXiv:1010.0202).
  //      Enlarge the ng x ng shape correlator with the Hankel offset ladder `offs` (=[0] for no
  //      Hankel): Big(t)[a,b] = matCA[t + offs[a] + offs[b]] (nb*ng x nb*ng).  Then metric = Big(t0fix)
  //      FIXED for all t; rebase at trebase fixes the interpolator V; reduced Cr(t)=V^T Big(t) V,
  //      whitened at t0fix, gives lambda_s(t); effmass = log(lam(t)/lam(t+1)) CONSECUTIVE (dt=1),
  //      dimensionless a_t*m when at=1.  Contrast t0mode=0 = SLIDING metric C(t) below (original).
  if(t0mode == 1){
    const int nb = (int)offs.size();
    const int nbig = nb*ng;
    int omax = 0;
    for(int o : offs) omax = std::max(omax, o);
    const int tmax = Nthalf - 2*omax;                  // Big(t) valid for t = 0 .. tmax-1
    for(int ti=0; ti<(int)out.size(); ti++)
      for(int s=0; s<nstates_g; s++) out[ti][s] = std::nan("");
    if(tmax < 2 || t0fix >= tmax || trebase >= tmax) return;   // not enough range -> all nan
    std::vector<MatrixXd> Big(tmax);
    for(int t=0; t<tmax; t++){
      MatrixXd B = MatrixXd::Zero(nbig, nbig);
      for(int a=0; a<nb; a++)
        for(int b=0; b<nb; b++)
          B.block(a*ng, b*ng, ng, ng) = matCA[t + offs[a] + offs[b]];
      Big[t] = 0.5*(B + B.transpose());
    }
    MatrixXd si0 = inv_sqrt_sym(Big[t0fix], rtol, kd);
    MatrixXd Creb = si0 * Big[trebase] * si0;          // rebase: Big(t0)^{-1/2} Big(trebase) Big(t0)^{-1/2}
    Eigen::SelfAdjointEigenSolver<MatrixXd> esr(0.5*(Creb+Creb.transpose()));
    MatrixXd Ur = esr.eigenvectors();                  // ascending
    const int nn = (int)esr.eigenvalues().size();
    MatrixXd Vop(nbig, nstates_g);
    for(int s=0; s<nstates_g; s++) Vop.col(s) = si0 * Ur.col(nn-1-s);  // top (descending) nstates_g
    MatrixXd Cr0 = Vop.transpose() * Big[t0fix] * Vop;
    MatrixXd sr0 = inv_sqrt_sym(Cr0, rtol, kd);
    std::vector<VectorXd> lam(tmax);
    for(int ti=1; ti<tmax; ti++){
      MatrixXd Mr = sr0 * (Vop.transpose()*Big[ti]*Vop) * sr0;
      Eigen::SelfAdjointEigenSolver<MatrixXd> es(0.5*(Mr+Mr.transpose()));
      lam[ti] = es.eigenvalues();                      // ascending
    }
    for(int ti=1; ti<tcut; ti++){
      if(ti+1 > tmax-1) break;
      const int nl = (int)lam[ti].size();
      for(int s=0; s<nstates_g; s++){
        const int idx = nl-1-s;                        // descending: s=0 = largest = ground
        out[ti-1][s] = std::log(lam[ti](idx) / lam[ti+1](idx)) / (dt*at);
      }
    }
    return;
  }

  // ---- SLIDING-t0 GEVP (t0mode=0, DEFAULT, ORIGINAL): metric = C(t) at each t ----
  MatrixXd Vre;
  if(use_rebase){
    MatrixXd si0 = inv_sqrt_sym(matCA[trebase], rtol, kd);
    MatrixXd M0 = si0 * matCA[trebase+dt] * si0;
    Eigen::SelfAdjointEigenSolver<MatrixXd> es0(M0);
    Vre = es0.eigenvectors();
  }
  for(int ti=1; ti<tcut; ti++){
    MatrixXd si = inv_sqrt_sym(matCA[ti], rtol, kd);
    MatrixXd M = si * matCA[ti+dt] * si;
    MatrixXd Mt;
    if(use_rebase){
      MatrixXd Mrot = Vre.transpose() * M * Vre;
      const int off = omit_top ? (ng-1-nstates_g) : (ng-nstates_g);
      Mt = Mrot.block(off, off, nstates_g, nstates_g);
    } else {
      Mt = M;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(Mt);
    VectorXd ev = es.eigenvalues();
    const int nev = (int)ev.size();
    for(int s=0; s<nstates_g; s++) out[ti-1][s] = -std::log(ev(nev - nstates_g + s)) / (dt*at);
  }
}

int main(int argc, char* argv[]){
  if(argc < 7){
    std::cerr << "usage: " << argv[0]
              << " <dir> <N_FLOW> <NCH> <drop_l0> <at> <outfile> [nops2=8] [tcut=5] [binsize=0auto] [kmin=100] [rtol=1e-8] [prefix=F_corr] [trebase=1]\n";
    return 1;
  }
  const std::string dir = argv[1];
  const int N_FLOW = std::atoi(argv[2]);
  const int NCH = std::atoi(argv[3]);
  const int drop_l0 = std::atoi(argv[4]);
  const double at = std::atof(argv[5]);
  const std::string outfile = argv[6];
  const int nops2 = (argc>7)  ? std::atoi(argv[7])  : 8;
  const int tcut  = (argc>8)  ? std::atoi(argv[8])  : 5;
  int binsize     = (argc>9)  ? std::atoi(argv[9])  : 0;
  const int kmin  = (argc>10) ? std::atoi(argv[10]) : 100;
  const double rtol = (argc>11) ? std::atof(argv[11]) : 1.0e-8; // GEVP-metric pruning threshold
  const std::string prefix = (argc>12) ? std::string(argv[12]) : std::string("F_corr"); // h5 file prefix
  const int trebase = (argc>13) ? std::atoi(argv[13]) : 1; // GEVP rebase/metric folded-time index t0
  // --- diagnostic knobs (default = original behavior) ---
  const std::string orbits_arg = (argc>14) ? std::string(argv[14]) : std::string(""); // iorbit subset, e.g. "0" or "0,3"; empty=all
  const int vacsub = (argc>15) ? std::atoi(argv[15]) : 1; // 1: subtract <O_i><O_j> (connected); 0: raw
  const int keepl0 = (argc>16) ? std::atoi(argv[16]) : 1; // 0: additionally drop all l=0 ops (even if present in h5)
  const double evec_t = (argc>17) ? std::atof(argv[17]) : 0.0; // >0: dump ground GEVP eigenvector (op basis) at this t
  const std::string lsel_arg = (argc>18) ? std::string(argv[18]) : std::string(""); // comma l's to KEEP (e.g. "1"); empty=all
  const int use_rebase = (argc>19) ? std::atoi(argv[19]) : 1; // 0: raw full GEVP (no t0 rebase / no truncation subspace)
  const int omit_top = (argc>20) ? std::atoi(argv[20]) : 0;   // 1: omit the top (largest-lambda) rebase mode (drop the l=0 constant/vacuum)
  const int deg_avg = (argc>21) ? std::atoi(argv[21]) : 0;    // degenerate l=1 m-triplet combination: 0=per-state output; 1=UNIT-weight mean; 2=INVERSE-VARIANCE weighted mean (weights from full-sample per-state jackknife variances, fixed, applied per sample then jackknifed -> correct error through the state correlations)
  const std::string jkdump = (argc>22) ? std::string(argv[22]) : std::string(""); // non-empty: dump the per-(t,state) jackknife-sample effmasses to this file for an external correlated constant fit
  const int t0mode = (argc>24) ? std::atoi(argv[24]) : 0;    // GEVP t0 handling: 0=SLIDING metric C(t) (default, original); 1=FIXED-t0 metric C(t0fix), effmass=log(lam_t/lam_{t+1})
  const int t0fix = (argc>25) ? std::atoi(argv[25]) : 0;     // fixed-t0 metric folded-sep INDEX (t0mode=1 only); rebase index = trebase (arg13)
  const std::string offs_arg = (argc>26) ? std::string(argv[26]) : std::string("0"); // Hankel offset ladder (t0mode=1), comma list e.g. "0,2"; "0"=no Hankel (shape basis only)
  std::vector<int> offs;
  {
    std::stringstream ss(offs_arg);
    std::string tok;
    while(std::getline(ss, tok, ',')){
      if(!tok.empty()) offs.push_back(std::atoi(tok.c_str()));
    }
    if(offs.empty()) offs.push_back(0);
  }
  const int per_m = (argc>23) ? std::atoi(argv[23]) : 0;     // 1: EXPLICIT per-m GEVP -- run a SEPARATE GEVP on each m-block (the n_shapes shapes of a fixed (l,m)); nops2 = states kept PER m-block; total reported states = (#m in the l-sector) * nops2 (e.g. F l=1: nops2=1 -> the 3 m-grounds). 0: one pooled GEVP over all (orbit,m) ops.

  const int Nt = 128;
  const int dt = 1;
  // Nthalf (number of folded time separations) is set AFTER the h5 is inspected: it depends on
  // whether the driver stored the FULL Nt separations (legacy density files) or only a small
  // window dt = 0..TMAX (shapes files). See the ntsep detection below.
  int Nthalf = Nt/2 - 1;   // folded dt = 1 .. Nt/2-1 (legacy default; overwritten for windowed h5)
  bool windowed = false;   // true when the h5 stores only dt = 0..ntsep-1 with ntsep < Nt

  // the (l,m) channel list is built AFTER inspecting the h5 (the driver writes n_lm, so the l-tower
  // auto-adapts, e.g. l=0..2 for the disk-reduced production or l=0..3 for older files).

  // ---- gather config indices k >= kmin ----
  std::vector<int> ks;
  for(const auto& e : std::filesystem::directory_iterator(dir)){
    const std::string fn = e.path().filename().string();
    const std::string pfx = prefix + ".";
    if(fn.rfind(pfx,0)==0 && fn.size()>3 && fn.substr(fn.size()-3)==".h5"){
      const int k = std::stoi(fn.substr(pfx.size())); // "<prefix>.123.h5" -> 123
      if(k>=kmin) ks.push_back(k);
    }
  }
  std::sort(ks.begin(), ks.end());
  const int Nc = ks.size();
  if(Nc==0){ std::cerr << "no F_corr files with k>=" << kmin << " in " << dir << "\n"; return 2; }

  // ---- nops + stored-separation count AUTO-DETECTED from the h5 (works for BOTH density and
  //      shapes layouts; op%n_lm=ilm). dims = {ntsep, nops*nops}. ntsep < Nt => windowed storage
  //      (driver kept only dt = 0..ntsep-1); the backward fold leg is then the transpose of a
  //      stored forward slice via C_ij(Nt-d) = C_ji(d). ntsep == Nt => legacy full-Nt storage. ----
  int nops = 0;
  int ntsep = Nt;
  int n_lm_full = 0;   // driver's measured (l,m) count; 0 for old h5 that predate the n_lm dataset
  int n_shapes_h5 = 0; // driver's shape/orbit count (block layout only)
  bool blocked = false;  // new drivers store only the per-(l,m) shape blocks ("F_corr_blk")
  std::vector<int> ell_h5, em_h5;   // explicit (l,m) list from the h5 (new drivers); empty for old h5
  {
    HighFive::File hf(dir+"/"+prefix+"."+std::to_string(ks[0])+".h5", HighFive::File::ReadOnly);
    if(hf.exist("n_lm")){ std::vector<int> tmp; hf.getDataSet("n_lm").read(tmp); n_lm_full = tmp[0]; }
    if(hf.exist("n_shapes")){ std::vector<int> tmp; hf.getDataSet("n_shapes").read(tmp); n_shapes_h5 = tmp[0]; }
    else if(hf.exist("n_orbits")){ std::vector<int> tmp; hf.getDataSet("n_orbits").read(tmp); n_shapes_h5 = tmp[0]; } // legacy key
    if(hf.exist("ell") && hf.exist("em")){ hf.getDataSet("ell").read(ell_h5); hf.getDataSet("em").read(em_h5); }
    if(hf.exist("F_corr_blk")){
      // per-(l,m) block layout: F_corr_blk = (nsep, n_lm*n_shapes^2); reconstruct nops = n_shapes*n_lm.
      blocked = true;
      const std::vector<size_t> dims = hf.getDataSet("F_corr_blk").getDimensions();
      ntsep = (int)dims.front();
      nops = n_shapes_h5 * n_lm_full;
      const long want = (long)n_lm_full * n_shapes_h5 * n_shapes_h5;
      if((long)dims.back() != want){ std::cerr << "F_corr_blk ncol " << dims.back() << " != n_lm*n_shapes^2="
                                               << want << "\n"; return 4; }
    } else {
      const std::vector<size_t> dims = hf.getDataSet("F_corr").getDimensions();
      ntsep = (int)dims.front();
      const long ncol0 = (long)dims.back();
      nops = (int)std::lround(std::sqrt((double)ncol0));
      if((long)nops*nops != ncol0){ std::cerr << "F_corr ncol " << ncol0 << " not a perfect square\n"; return 4; }
    }
  }
  const int ncol = nops*nops;

  // ---- (l,m) list per ilm channel ----
  // Preferred: read the EXACT list from the h5 ell/em datasets (new drivers; supports non-contiguous
  // l-towers, e.g. F saves l=1,2 and F^2 saves l=0,1). Legacy fallback: reconstruct l=0..lmax from n_lm.
  std::vector<std::pair<int,int>> lm;
  if(!ell_h5.empty() && ell_h5.size()==em_h5.size()){
    for(size_t i=0; i<ell_h5.size(); i++) lm.push_back({ell_h5[i], em_h5[i]});
  } else {
    int lmax = 3;
    if(n_lm_full > 0) lmax = (int)std::lround(std::sqrt((double)(drop_l0 ? n_lm_full+1 : n_lm_full))) - 1;
    for(int l=0; l<=lmax; l++){
      for(int m=-l; m<=l; m++){
        if(drop_l0 && l==0) continue;
        lm.push_back({l,m});
      }
    }
  }
  const int n_lm = (int)lm.size();
  if(nops % n_lm != 0){ std::cerr << "nops " << nops << " not a multiple of n_lm=" << n_lm
                                  << " (h5 n_lm=" << n_lm_full << ")\n"; return 4; }
  windowed = (ntsep < Nt);
  if(windowed) Nthalf = ntsep - 2;   // folds idx 0..ntsep-3 use forward slice tt and transpose of tt+1
  if(Nthalf < tcut){ std::cerr << "stored window too small: Nthalf=" << Nthalf << " < tcut=" << tcut
                               << " (need a larger TMAX_CORR in the driver)\n"; return 4; }
  if(trebase < 0 || trebase+dt > Nthalf-1){
    std::cerr << "trebase=" << trebase << " out of range [0," << (Nthalf-1-dt) << "]\n"; return 4;
  }
  if(t0mode==1 && (t0fix < 0 || t0fix > Nthalf-1)){
    std::cerr << "t0fix=" << t0fix << " out of range [0," << (Nthalf-1) << "]\n"; return 4;
  }
  (void)N_FLOW; (void)NCH; // superseded by h5 auto-detect; kept for CLI back-compat

  // opset (l,m) per operator: op mod n_lm = ilm
  // (density: iflow*(NCH*n_lm)+ich*n_lm+ilm ; shapes: iorbit*n_lm+ilm)
  std::vector<int> opL(nops), opM(nops);
  for(int op=0; op<nops; op++){
    opL[op] = lm[op % n_lm].first;
    opM[op] = lm[op % n_lm].second;
  }

  // operator subset kept for the GEVP: chosen orbits (shape decomposition) and, if keepl0==0, drop l=0.
  const int n_shapes = nops / n_lm;
  std::vector<bool> orbit_keep(n_shapes, orbits_arg.empty());
  if(!orbits_arg.empty()){
    std::stringstream ss(orbits_arg);
    std::string tok;
    while(std::getline(ss, tok, ',')){
      if(tok.empty()) continue;
      const int io = std::stoi(tok);
      if(io>=0 && io<n_shapes) orbit_keep[io] = true;
    }
  }
  // l-sector selection (keep the full l=0..3 tower in the DRIVER; focus one sector here). Empty = all l.
  std::vector<bool> l_keep(4, lsel_arg.empty());
  if(!lsel_arg.empty()){
    std::stringstream ls(lsel_arg);
    std::string tok;
    while(std::getline(ls, tok, ',')){
      if(tok.empty()) continue;
      const int l = std::stoi(tok);
      if(l>=0 && l<=3) l_keep[l] = true;
    }
  }
  std::vector<int> keep;
  for(int op=0; op<nops; op++){
    if(!orbit_keep[op / n_lm]) continue;
    if(!keepl0 && opL[op]==0) continue;
    if(!l_keep[opL[op]]) continue;
    keep.push_back(op);
  }
  const int nk = (int)keep.size();
  if(nk==0){ std::cerr << "operator subset empty\n"; return 4; }
  std::cout << "# subset: n_shapes=" << n_shapes << " kept_ops=" << nk
            << " vacsub=" << vacsub << " keepl0=" << keepl0 << " orbits=[" << orbits_arg
            << "] lsel=[" << lsel_arg << "]\n";

  // GEVP grouping. per_m=1: partition the kept ops into one m-block per (l,m) (each = the n_shapes
  // shapes at that m) and run a SEPARATE GEVP on each -- correct for the exactly-degenerate l=1,2
  // multiplets (blocks are decoupled by symmetry). per_m=0: one pooled GEVP over all kept ops.
  std::vector<std::vector<int>> groups;
  if(per_m){
    std::vector<int> mvals;
    for(int a=0; a<nk; a++){
      const int m = opM[keep[a]];
      if(std::find(mvals.begin(), mvals.end(), m) == mvals.end()) mvals.push_back(m);
    }
    std::sort(mvals.begin(), mvals.end());
    for(size_t g=0; g<mvals.size(); g++){
      std::vector<int> grp;
      for(int a=0; a<nk; a++){ if(opM[keep[a]]==mvals[g]) grp.push_back(keep[a]); }
      groups.push_back(grp);
    }
  } else {
    groups.push_back(keep);
  }
  const int ngroups = (int)groups.size();
  const int gsize = (int)groups[0].size();
  std::cout << "# GEVP groups: ngroups=" << ngroups << " gsize=" << gsize
            << " (per_m=" << per_m << ")\n";

  if(binsize<=0) binsize = std::max(1, Nc/50);
  const int nbins = Nc/binsize;
  const int N = nbins*binsize; // configs used (drop remainder), as in the notebook Jackknife

  std::cout << "# dir=" << dir << " nops=" << nops << " Nc=" << Nc
            << " binsize=" << binsize << " nbins=" << nbins << " N=" << N << std::endl;

  // ---- accumulators: total + per-bin sums of folded corr (ncol x Nthalf) and one-point (nops) ----
  std::vector<MatrixXd> Sb(nbins, MatrixXd::Zero(ncol, Nthalf));
  MatrixXd Stot = MatrixXd::Zero(ncol, Nthalf);
  std::vector<VectorXd> Sb1(nbins, VectorXd::Zero(nops));
  VectorXd Stot1 = VectorXd::Zero(nops);

  std::vector<double> row(ncol);
  for(int ci=0; ci<N; ci++){
    const int k = ks[ci];
    const int b = ci/binsize;

    // read F_corr (Nt x ncol) and one-point F (nops) from the per-config h5
    HighFive::File hf(dir+"/"+prefix+"."+std::to_string(k)+".h5", HighFive::File::ReadOnly);
    std::vector<std::vector<double>> Fc;
    if(blocked){
      // read the per-(l,m) blocks and EXPAND to the block-diagonal full matrix (cross-(l,m) = 0).
      // blk[tt][ilm*nob^2 + a*nob + b] -> Fc[tt][(a*n_lm+ilm)*nops + (b*n_lm+ilm)].
      std::vector<std::vector<double>> blk;
      hf.getDataSet("F_corr_blk").read(blk);
      Fc.assign(ntsep, std::vector<double>(ncol, 0.0));
      for(int tt=0; tt<ntsep; tt++){
        for(int ilm=0; ilm<n_lm; ilm++){
          for(int a=0; a<n_shapes; a++){
            const int oi = a*n_lm + ilm;
            for(int b=0; b<n_shapes; b++){
              const int oj = b*n_lm + ilm;
              Fc[tt][oi*nops+oj] = blk[tt][ilm*n_shapes*n_shapes + a*n_shapes + b];
            }
          }
        }
      }
    } else {
      hf.getDataSet("F_corr").read(Fc);
    }
    std::vector<double> f1v;
    hf.getDataSet("F").read(f1v);

    // fold pos/neg t: fold[tt-1] = C(tt) + C(Nt-1-tt), tt=1..Nthalf. The backward leg C(Nt-1-tt)
    // = C_ji(tt+1) = transpose of the forward slice Fc[tt+1] (periodicity C_ij(Nt-d)=C_ji(d)).
    //  - legacy full-Nt h5: read the physical backward slice Fc[Nt-1-tt] directly (bit-identical).
    //  - windowed h5 (only dt=0..ntsep-1 stored): reconstruct it as Fc[tt+1]^T.
    for(int tt=1; tt<=Nthalf; tt++){
      const int idx = tt-1;
      const double* rup = Fc[tt].data();
      if(windowed){
        const double* rtp = Fc[tt+1].data();   // backward leg = transpose of stored forward slice
        for(int i=0; i<nops; i++){
          for(int j=0; j<nops; j++){
            const double v = rup[i*nops+j] + rtp[j*nops+i];
            Sb[b](i*nops+j,idx) += v;
            Stot(i*nops+j,idx) += v;
          }
        }
      } else {
        const double* rdn = Fc[Nt-1-tt].data();
        for(int c=0; c<ncol; c++){
          const double v = rup[c] + rdn[c];
          Sb[b](c,idx) += v;
          Stot(c,idx) += v;
        }
      }
    }
    for(int op=0; op<nops; op++){
      Sb1[b](op) += f1v[op];
      Stot1(op) += f1v[op];
    }

    if((ci%500)==0) std::cout << "# read " << ci << "/" << N << std::endl;
  }

  // ---- jackknife over bins ----
  // states kept per group (nops2, capped by group size), and total across groups.
  const int nstates_g = omit_top ? std::min(nops2, gsize-1) : std::min(nops2, gsize);
  const int nstates = ngroups * nstates_g;
  const int ntpts = tcut-1;
  // jk[b][ti][state]
  std::vector<std::vector<std::vector<double>>> jk(
      nbins, std::vector<std::vector<double>>(ntpts, std::vector<double>(nstates, 0.0)));

  // ---- UNIFORM metric pruning (NM 2026-10-12): decide every kept dimension ONCE on the FULL sample ----
  std::vector<KeepDims> kdims(ngroups);
  {
    MatrixXd conn_full = Stot / (double)N;
    VectorXd m1_full = Stot1 / (double)N;
    prepare_conn(conn_full, m1_full, nops, Nthalf, vacsub, opL, opM);
    for(int g=0; g<ngroups; g++){
      std::vector<std::vector<double>> gout(ntpts, std::vector<double>(nstates_g, 0.0));
      kdims[g].determine = true;
      gevp_group(groups[g], conn_full, nops, Nthalf, trebase, dt, tcut, rtol, use_rebase, omit_top, at, nstates_g, gout, t0mode, t0fix, offs, kdims[g]);
      kdims[g].determine = false;
      std::cout << "# uniform pruning: group " << g << " (size " << groups[g].size() << ") kept dims (full sample, rtol "
                << rtol << ") =";
      for(int nk : kdims[g].n) std::cout << " " << nk;
      std::cout << std::endl;
    }
  }

  for(int b=0; b<nbins; b++){
    const double denom = N - binsize;
    MatrixXd conn = (Stot - Sb[b]) / denom;      // (ncol x Nthalf)
    VectorXd m1   = (Stot1 - Sb1[b]) / denom;    // (nops)

    // vacuum subtraction C_ij -= <O_i><O_j> (vacsub) + zero_noise of rotationally forbidden entries
    prepare_conn(conn, m1, nops, Nthalf, vacsub, opL, opM);

    // per-group GEVP: run a separate GEVP on each group (one m-block when per_m, else the whole kept
    // subset), keeping nstates_g states each; write into jk at slots [g*nstates_g + s]. jk keeps the
    // RAW per-state effmasses; the degenerate-multiplet combination (unit / inverse-variance) is
    // formed in the write section below so the weights can use the full-sample per-state variances.
    // The metric pruning REPLAYS the full-sample kept dimensions (kdims[g]) identically in every sample.
    for(int g=0; g<ngroups; g++){
      std::vector<std::vector<double>> gout(ntpts, std::vector<double>(nstates_g, 0.0));
      kdims[g].rewind();
      gevp_group(groups[g], conn, nops, Nthalf, trebase, dt, tcut, rtol, use_rebase, omit_top, at, nstates_g, gout, t0mode, t0fix, offs, kdims[g]);
      for(int ti=0; ti<ntpts; ti++){
        for(int s=0; s<nstates_g; s++) jk[b][ti][g*nstates_g + s] = gout[ti][s];
      }
    }
  }
  for(int g=0; g<ngroups; g++){
    std::cout << "# uniform pruning check: group " << g << ": per-sample rtol rule would differ in "
              << kdims[g].n_mismatch << " (sample,call) pairs; non-positive kept eigenvalues (-> NaN samples, LOUD): "
              << kdims[g].n_nonpos << std::endl;
  }

  // ---- optional: dump per-(t,state) jackknife samples for an external correlated constant fit ----
  // Layout: line 1 "# nbins ntpts nstates at"; then nbins rows, each with ntpts*nstates values,
  // column index = ti*nstates + s (ti=0..ntpts-1, s=0..nstates-1). The fitter builds the (ntpts*nstates)
  // covariance from these bins (degenerate l=1: all 3 states share one mass -> fit them simultaneously).
  if(!jkdump.empty()){
    std::ofstream jf(jkdump);
    jf << std::scientific << std::setprecision(12);
    jf << "# " << nbins << " " << ntpts << " " << nstates << " " << at << "\n";
    for(int b=0; b<nbins; b++){
      for(int ti=0; ti<ntpts; ti++){
        for(int s=0; s<nstates; s++) jf << jk[b][ti][s] << (ti==ntpts-1 && s==nstates-1 ? '\n' : ' ');
      }
    }
    jf.close();
    std::cout << "# wrote jackknife-sample dump " << jkdump << " (" << nbins << " bins x " << ntpts*nstates << " points)\n";
  }

  // ---- optional: dump the ground-state GEVP eigenvector (op basis) at t=evec_t ----
  // v solves C(t+dt) v = lambda C(t) v ; v_i = weight of operator O_i in the optimal interpolator.
  // Reported sorted by |v_i| with (orbit/shape, l, m) labels so one can see which ops the state overlaps.
  if(evec_t > 0.0){
    const int ti = (int)std::lround(evec_t/at);
    if(ti < 1 || ti+dt > Nthalf-1){
      std::cerr << "evec_t=" << evec_t << " -> ti=" << ti << " out of range [1," << (Nthalf-1-dt) << "]\n";
    } else {
      MatrixXd conn = Stot / (double)N;   // full-sample connected correlator (ncol x Nthalf)
      VectorXd m1 = Stot1 / (double)N;
      if(vacsub){
        for(int i=0; i<nops; i++){
          for(int j=0; j<nops; j++){
            const double p = m1(i)*m1(j);
            for(int t=0; t<Nthalf; t++) conn(i*nops+j,t) -= p;
          }
        }
      }
      for(int i=0; i<nops; i++){
        for(int j=0; j<nops; j++){
          const int l1=opL[i], l2=opL[j], mm1=opM[i], mm2=opM[j];
          bool z=false;
          if((l1==1 && l2==0) || (l1==0 && l2==1)) z=true;
          else if(l1==1 && l2==1 && mm1!=mm2) z=true;
          if(z) for(int t=0; t<Nthalf; t++) conn(i*nops+j,t) = 0.0;
        }
      }
      MatrixXd Ca(nk,nk), Cb(nk,nk);
      for(int a=0; a<nk; a++){
        for(int b=0; b<nk; b++){
          Ca(a,b) = conn(keep[a]*nops+keep[b], ti);
          Cb(a,b) = conn(keep[a]*nops+keep[b], ti+dt);
        }
      }
      Ca = 0.5*(Ca + Ca.transpose());
      Cb = 0.5*(Cb + Cb.transpose());
      KeepDims kd_evec;
      MatrixXd si = inv_sqrt_sym(Ca, rtol, kd_evec);
      MatrixXd M = si * Cb * si;
      M = 0.5*(M + M.transpose());
      Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
      VectorXd ev = es.eigenvalues();     // ascending -> ground = largest = last
      const char* SHAPES[7] = {"triangle","rect","fig8","three-tri","star","trio","five-six"}; // consolidated 7-shape basis (all L)
      std::cout << std::fixed << std::setprecision(4);
      // dump the top `nstates` GEVP modes (rev=0 = ground/largest-lambda; for the non-subtracted l=0
      // F^2 channel rev=0 is the constant/vacuum mode Delta~0 and rev=1 is the physical 0++).
      const int nev_dump = std::min(nstates, nk);
      for(int rev=0; rev<nev_dump; rev++){
        const int col = nk-1-rev;          // ground first, then next-lightest, ...
        VectorXd w = es.eigenvectors().col(col);
        VectorXd v = si * w;               // op-basis generalized eigenvector
        v.normalize();
        int imax=0;
        for(int a=1; a<nk; a++){ if(std::abs(v(a)) > std::abs(v(imax))) imax=a; }
        if(v(imax) < 0.0) v = -v;          // fix overall sign: dominant component positive
        // physical overlap Z_i = <0|O_i|n> propto (C(t0) v)_i (covariant partner of the interpolator
        // weight v). Normalized to unit norm so |Z| is directly comparable to |v|.
        VectorXd Zv = Ca * v;
        Zv.normalize();
        if(Zv(imax) < 0.0) Zv = -Zv;
        std::vector<std::pair<double,int>> mag(nk);
        for(int a=0; a<nk; a++) mag[a] = std::make_pair(std::abs(v(a)), a);
        std::sort(mag.rbegin(), mag.rend()); // descending by |v_i|
        std::cout << "# eigenvector mode " << rev << " at t=" << (ti*at)
                  << "  Delta=" << (-std::log(ev(col))/(dt*at))
                  << "  lambda=" << ev(col) << "  nk=" << nk << "\n";
        std::cout << "# rank    |v|      v        |Z|      Z      orbit  shape       l   m   (Z=C(t0)v, both unit-norm)\n";
        for(int r=0; r<nk; r++){
          const int a = mag[r].second;
          const int op = keep[a];
          const int iorbit = op / n_lm;
          const int ilm = op % n_lm;
          const char* sh = (n_shapes==7 && iorbit<7) ? SHAPES[iorbit] : "orbit"; // consolidated 7-shape basis
          std::cout << "  " << std::setw(3) << r
                    << "  " << std::setw(8) << std::abs(v(a))
                    << "  " << std::setw(8) << std::showpos << v(a) << std::noshowpos
                    << "  " << std::setw(8) << std::abs(Zv(a))
                    << "  " << std::setw(8) << std::showpos << Zv(a) << std::noshowpos
                    << "   " << std::setw(3) << iorbit
                    << "  " << std::setw(9) << sh
                    << "  " << std::setw(2) << lm[ilm].first
                    << "  " << std::setw(2) << lm[ilm].second << "\n";
        }
      }
      std::cout << std::defaultfloat;
    }
  }

  // ---- jackknife mean / error ----
  // out row: t_phys  ground_mean ground_err  [state s: mean err ...]  (ground = state nstates-1)
  // deg_avg 1/2: report ONLY the triplet combination (single reported series); ground = the average.
  std::ofstream out(outfile);
  out << std::scientific << std::setprecision(10);
  if(deg_avg==1)
    out << "# t   avg_mean avg_err   (UNIT-weight mean of " << nstates << " l=1 effmasses per jk sample; ground=avg)\n";
  else if(deg_avg==2)
    out << "# t   wavg_mean wavg_err   (INVERSE-VARIANCE weighted mean of " << nstates << " l=1 effmasses; fixed weights from full-sample per-state variances; ground=wavg)\n";
  else
    out << "# t   ground_mean ground_err   then per-state (s=0..nstates-1) mean err ; ground=s" << (nstates-1) << "\n";
  std::cout << std::scientific << std::setprecision(6);
  std::cout << "# ground-state Delta_eff(t):" << std::endl;
  for(int ti=0; ti<ntpts; ti++){
    const double tphys = (ti+1)*at;
    // per-state jackknife mean/err (always computed; used directly for deg_avg=0 and as
    // the inverse-variance weights for deg_avg=2)
    std::vector<double> mean(nstates,0.0), err(nstates,0.0);
    for(int s=0; s<nstates; s++){
      double mu=0.0;
      for(int b=0; b<nbins; b++) mu += jk[b][ti][s];
      mu /= nbins;
      double var=0.0;
      for(int b=0; b<nbins; b++) var += (jk[b][ti][s]-mu)*(jk[b][ti][s]-mu);
      var *= (double)(nbins-1)/nbins;
      mean[s]=mu;
      err[s]=std::sqrt(var);
    }

    if(deg_avg){
      // fixed weights w_s: unit (deg_avg=1) or inverse-variance (deg_avg=2), normalized to sum 1
      std::vector<double> w(nstates, 1.0/nstates);
      if(deg_avg==2){
        double wsum=0.0;
        for(int s=0; s<nstates; s++){
          w[s] = (err[s]>0.0) ? 1.0/(err[s]*err[s]) : 0.0;
          wsum += w[s];
        }
        for(int s=0; s<nstates; s++) w[s] /= wsum;
      }
      // apply the fixed weights per jackknife sample, then jackknife the linear combination
      // (captures the state-state correlations in the error, unlike 1/sqrt(sum 1/var))
      double Mmu=0.0;
      std::vector<double> Mb(nbins,0.0);
      for(int b=0; b<nbins; b++){
        double m=0.0;
        for(int s=0; s<nstates; s++) m += w[s]*jk[b][ti][s];
        Mb[b]=m;
        Mmu+=m;
      }
      Mmu/=nbins;
      double Mvar=0.0;
      for(int b=0; b<nbins; b++) Mvar += (Mb[b]-Mmu)*(Mb[b]-Mmu);
      Mvar *= (double)(nbins-1)/nbins;
      const double Merr=std::sqrt(Mvar);
      out << tphys << " " << Mmu << " " << Merr << " " << Mmu << " " << Merr << "\n";
      std::cout << "#   t=" << tphys << "  m_ground=" << Mmu << " +/- " << Merr << std::endl;
    } else {
      const int g = nstates-1;
      out << tphys << " " << mean[g] << " " << err[g];
      for(int s=0; s<nstates; s++) out << " " << mean[s] << " " << err[s];
      out << "\n";
      std::cout << "#   t=" << tphys << "  m_ground=" << mean[g] << " +/- " << err[g] << std::endl;
    }
  }
  out.close();
  std::cout << "# wrote " << outfile << std::endl;
  return 0;
}
