// jj_commut_tp_claude.cc
// -----------------------------------------------------------------------------
// Temporal (tp) conserved-current correlator from the COMMUTATOR shortcut, using the dense D_ov and
// P = D_m^{-1} already produced by jj_propagator_deter.  NO overlap multishift, NO matmul -- a single
// O(N^2) pass.  (jj_exact_freefield_impl_plan_claude.md, "COMMUTATOR shortcut".)
//
// Identity: the conserved current contracted with a gauge twist Theta=diag(theta_x) is exact,
//   [D_ov, Theta] = sum_{wz} (theta_w - theta_z) k^{wz}        (gauge covariance; check_conserved_current).
// The temporal projection is the time-step theta_x = 1[t_x > t]: on the time circle this is a TWO-edge
// cut (at t and at the wrap); since D_ov (hence K) is local in t over ~10-30 slices and the wrap edge is
// ~Nt/2 away, the second edge is suppressed by ~1e-8 at Nt=128, so K^{tp}(t) = [D_ov, Pi_{>t}] is a clean
// single cut.  Massless (D_ov P = I), with A(t)=K^{tp}(t)P = G_t - Pi_{>t}, G_t=D_ov Pi_{>t} D_ov^{-1}:
//   disc(t) = tr A(t) = 0   (exactly; computed here as a PDov-diagonal sanity check),
//   conn(t0,t) = tr[A(t0)A(t)] = 2 tr(Pi_{>max(t0,t)}) - S(t0,t) - S(t,t0),
//   S(a,b) = sum_{t_i>a, t_j>b} (D_ov)_{ij} P_{ji}     (masked sum of M_{ij}=(D_ov)_{ij}P_{ji}).
// We block-sum M over the Nt time-blocks (m[u][v]) once, take its 2D suffix sum (S), and every conn is
// O(1).  conn(t0,t) = C(t0-t) + (Q^2-type constant) + ~1e-8 wrap terms.
//
// CLI: --ens-dir(omit=free) --mass-re --mass-im --L --n-t0 --ninter.
//   reads data_<ESNID>/prop_deter_L<L>/Dinv.<k>.h5 -> writes data_<ESNID>/corr_deter_commut_L<L>/corr.<k>.h5
// Build: g++ -std=c++17 (HighFive header-only + hdf5).  Pure host; no CUDA.
// -----------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <getopt.h>
#include <highfive/H5File.hpp>

using Idx = std::int64_t;
using Complex = std::complex<double>;

static std::string basename_of(std::string s){
  if(!s.empty() && s.back()=='/') s.pop_back();
  auto p = s.find_last_of('/');
  if(p!=std::string::npos) s = s.substr(p+1);
  return s;
}

static int read_scalar_int(HighFive::File& f, const std::string& key){
  std::vector<int> v;
  f.getDataSet(key).read(v);
  return v.at(0);
}

static void load_cmat(HighFive::File& f, const std::string& key, std::vector<Complex>& M){
  std::vector<double> re, im;
  f.getDataSet(key+"/real").read(re);
  f.getDataSet(key+"/imag").read(im);
  M.resize(re.size());
  for(size_t i=0;i<re.size();i++) M[i] = Complex(re[i], im[i]);
}

static void write_cvec(HighFive::File& f, const std::string& key, const std::vector<Complex>& C){
  std::vector<double> re(C.size()), im(C.size());
  for(size_t i=0;i<C.size();i++){ re[i]=C[i].real(); im[i]=C[i].imag(); }
  f.createDataSet(key+"/real", re);
  f.createDataSet(key+"/imag", im);
}

struct Args{ std::string ens_dir; double mass_re=0.0, mass_im=0.0; int L=1, n_t0=2, ninter=10; };

static void PrintHelp(){
  std::printf("jj_commut_tp: --ens-dir(omit=free) --mass-re --mass-im --L --n-t0 --ninter\n");
}

static void ParseArgs(int argc, char** argv, Args& a){
  static struct option lo[] = {
    {"ens-dir",required_argument,0,'e'}, {"mass-re",required_argument,0,'r'},
    {"mass-im",required_argument,0,'i'}, {"L",required_argument,0,'L'},
    {"n-t0",required_argument,0,'T'},    {"ninter",required_argument,0,'I'},
    {"help",no_argument,0,'h'}, {0,0,0,0} };
  int opt, idx;
  while((opt=getopt_long(argc,argv,"e:r:i:L:T:I:h",lo,&idx))!=-1){
    switch(opt){
      case 'e': a.ens_dir=optarg; break;
      case 'r': a.mass_re=std::stod(optarg); break;
      case 'i': a.mass_im=std::stod(optarg); break;
      case 'L': a.L=std::stoi(optarg); break;
      case 'T': a.n_t0=std::stoi(optarg); break;
      case 'I': a.ninter=std::stoi(optarg); break;
      case 'h': default: PrintHelp(); std::exit(0);
    }
  }
}

// block-sum M_{ij}=(D_ov)_{ij}P_{ji} over the Nt time-blocks: m[u*Nt+v] = sum_{t_i=u, t_j=v} M_{ij}.
// (D_ov, P are row-major length N^2; the time index of a flat index is index/Nx.)
static void block_sum_M(const std::vector<Complex>& Dov, const std::vector<Complex>& P,
                        Idx N, int Nt, Idx Nx, std::vector<Complex>& m){
  m.assign((size_t)Nt*Nt, Complex(0,0));
  for(Idx i=0;i<N;i++){
    const int u = (int)(i/Nx);
    const size_t row = (size_t)i*N;
    for(Idx j=0;j<N;j++){
      const int v = (int)(j/Nx);
      const Complex Mij = Dov[row+j] * P[(size_t)j*N+i];
      m[(size_t)u*Nt+v] += Mij;
    }
  }
}

// 2D suffix sum: Ssuf[a*(Nt+1)+b] = sum_{u>=a, v>=b} m[u][v]  (a,b in 0..Nt, the Nt row/col = 0).
static void suffix_sum(const std::vector<Complex>& m, int Nt, std::vector<Complex>& Ssuf){
  const int W = Nt+1;
  Ssuf.assign((size_t)W*W, Complex(0,0));
  for(int a=Nt-1;a>=0;a--){
    for(int b=Nt-1;b>=0;b--){
      Ssuf[(size_t)a*W+b] = m[(size_t)a*Nt+b]
        + Ssuf[(size_t)(a+1)*W+b] + Ssuf[(size_t)a*W+(b+1)] - Ssuf[(size_t)(a+1)*W+(b+1)];
    }
  }
}

int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  Args a; ParseArgs(argc, argv, a);
  const bool free_field = a.ens_dir.empty();
  const bool parity = (a.mass_im != 0.0);

  const std::string tag = free_field ? std::string("free") : basename_of(a.ens_dir);
  const std::string esnid = tag + "_vmRe" + std::to_string(a.mass_re) + "vmIm" + std::to_string(a.mass_im);
  const std::string propdir = "data_" + esnid + "/prop_deter_L" + std::to_string(a.L) + "/";
  const std::string outdir  = "data_" + esnid + "/corr_deter_commut_L" + std::to_string(a.L) + "/";
  std::filesystem::create_directories(outdir);
  std::cout << "# commutator-tp: esnid=" << esnid << (free_field?"  [free]":"  [interacting]") << "\n";

  const int k_lo = 0, k_step = free_field ? 1 : a.ninter, k_hi = free_field ? 1 : 1000000;
  for(int k=k_lo; k<k_hi; k+=k_step){
    const std::string pfile = propdir + "Dinv." + std::to_string(k) + ".h5";
    if(!std::filesystem::exists(pfile)){
      std::cout << "# k=" << k << " no propagator " << pfile << "\n";
      if(free_field) break; else continue;
    }
    const std::string h5path = outdir + "corr." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool c=false;
      try{ HighFive::File f(h5path, HighFive::File::ReadOnly); c=f.exist("complete"); }catch(...){}
      if(c){ std::cout << "# skip k=" << k << " (complete)\n"; if(free_field) break; else continue; }
    }

    Idx N; int Nt, n_sites; std::vector<Complex> Dov, P;
    {
      HighFive::File f(pfile, HighFive::File::ReadOnly);
      N = read_scalar_int(f, "N");
      Nt = read_scalar_int(f, "Nt");
      n_sites = read_scalar_int(f, "n_sites");
      if(!f.exist("Dov")){ std::cout << "# k=" << k << " no Dov in " << pfile << " (rerun propagator)\n"; if(free_field) break; else continue; }
      load_cmat(f, "Dov", Dov);
      load_cmat(f, "Dm_inv", P);
    }
    const Idx Nx = N / Nt;
    std::cout << "# k=" << k << "  N=" << N << " Nt=" << Nt << " Nx=" << Nx << "  loaded Dov,P\n";

    // M block-sum and 2D suffix sum
    std::vector<Complex> m, Ssuf;
    block_sum_M(Dov, P, N, Nt, Nx, m);
    suffix_sum(m, Nt, Ssuf);
    const int W = Nt+1;
    // S(t0,t) = Ssuf[(t0+1)*W + (t+1)] ;  tr(Pi_{>s}) = (Nt-1-s)*Nx
    // conn(t0,t) = 2*trPi(max) - S(t0,t) - S(t,t0)

    // disc(t) sanity = tr(G_t) - tr(Pi_{>t}) ; tr(G_t)=sum_{t_j>t}(P D_ov)_{jj}.  (=0 analytically.)
    std::vector<Complex> dj(N, Complex(0,0));   // d_j = (P D_ov)_{jj} = sum_k P_{jk}(D_ov)_{kj}
    for(Idx j=0;j<N;j++){
      const size_t rj = (size_t)j*N;
      Complex s(0,0);
      for(Idx kk=0;kk<N;kk++) s += P[rj+kk] * Dov[(size_t)kk*N+j];
      dj[j] = s;
    }
    std::vector<Complex> discvec(Nt, Complex(0,0));
    double max_disc = 0.0;
    for(int t=0;t<Nt;t++){
      Complex g(0,0);
      for(Idx j=0;j<N;j++){ if((int)(j/Nx) > t) g += dj[j]; }
      const double trPi = (double)(Nt-1-t) * (double)Nx;
      discvec[t] = g - Complex(trPi, 0.0);
      max_disc = std::max(max_disc, std::abs(discvec[t]));
    }

    // per-origin Vpp[b][dt] = conn(t0_b, (t0_b+dt)%Nt) ; and translation-averaged Gt_avg(dt)
    const int n_t0 = a.n_t0;
    std::vector<int> t0s(n_t0);
    for(int b=0;b<n_t0;b++) t0s[b] = b * (Nt/n_t0);
    std::vector<std::vector<Complex>> Vpp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
    std::vector<Complex> Gt_avg(Nt, Complex(0,0));
    for(int t0=0;t0<Nt;t0++){
      for(int t=0;t<Nt;t++){
        const int mx = (t0>t)?t0:t;
        const double trPi = (double)(Nt-1-mx) * (double)Nx;
        const Complex conn = Complex(2.0*trPi,0.0)
          - Ssuf[(size_t)(t0+1)*W + (t+1)] - Ssuf[(size_t)(t+1)*W + (t0+1)];
        const int dt = ((t - t0) % Nt + Nt) % Nt;
        Gt_avg[dt] += conn / (double)Nt;
      }
    }
    for(int b=0;b<n_t0;b++){
      const int t0 = t0s[b];
      for(int dt=0;dt<Nt;dt++){
        const int t = (t0+dt) % Nt;
        const int mx = (t0>t)?t0:t;
        const double trPi = (double)(Nt-1-mx) * (double)Nx;
        Vpp[b][dt] = Complex(2.0*trPi,0.0)
          - Ssuf[(size_t)(t0+1)*W + (t+1)] - Ssuf[(size_t)(t+1)*W + (t0+1)];
      }
    }

    // ---- atomic write data_<ESNID>/corr_deter_commut_L<L>/corr.<k>.h5 ----
    const std::string h5tmp = h5path + ".tmp";
    auto h5p = std::make_unique<HighFive::File>(h5tmp,
                 HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5 = *h5p;
    h5.createDataSet("N", std::vector<int>{(int)N});
    h5.createDataSet("Nt", std::vector<int>{Nt});
    h5.createDataSet("n_sites", std::vector<int>{n_sites});
    h5.createDataSet("t0s", t0s);
    h5.createDataSet("n_t0", std::vector<int>{n_t0});
    h5.createDataSet("parity", std::vector<int>{parity?1:0});
    write_cvec(h5, "h0/tp/Gt_avg", Gt_avg);
    for(int b=0;b<n_t0;b++) write_cvec(h5, "h0/t0_"+std::to_string(b)+"/tp/Vpp", Vpp[b]);
    write_cvec(h5, "h0/disc/tp/J", discvec);
    h5.createDataSet("complete", std::vector<int>{1});
    h5p.reset();
    std::filesystem::rename(h5tmp, h5path);

    std::cout << "#   wrote " << h5path << "   max|disc(=PDov-I diag)|=" << max_disc << "\n";
    std::cout << "#   Gt_avg(dt) dt=0,1,2,4,8,16,32,64: ";
    for(int dt : {0,1,2,4,8,16,32,64}) std::cout << Gt_avg[dt].real() << " ";
    std::cout << "\n";
    if(free_field) break;
  }
  return 0;
}
