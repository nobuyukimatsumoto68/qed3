// jj_contract_deter_claude.cu
// -----------------------------------------------------------------------------
// STAGE 2 (final) of the exact current-correlator validation: PURE LINEAR ALGEBRA.
// Loads the dense propagator P = D_m^{-1} (jj_propagator_deter, mass-dependent) and the dense
// projected current K^P(t) (jj_kbuild_exact, mass-INDEPENDENT) and contracts:
//   A^P  := K^P(0) . P            (cuBLAS Zgemm)
//   disc J^P(t) = tr(A^P(t)),   conn C^P(t0,t) = tr(A^P(t0) A^P(t)).
// FREE field (time-translation invariant): A^P(t)=S_t A^P(0) S_{-t}, so conn(dt)=tr(A0 . S_dt A0 S_-dt)
//   via the time-shift (layout idx = Nx*t + x, Nx=NS*N_SITES => time is the OUTER index).
// INTERACTING: no shift -> K^P(t) for all t is loaded (kbuild saved them) and conn=tr(A(t0) A(t)).
// NO overlap/Dirac/geometry here -- just cuBLAS + HDF5 (fast, re-runnable).
//
// L (= N_REFINE) COMPILE-TIME: -DN_REFINE_CLI=<L>.  CLI: --ens-dir (omit => free), --mass-re/--mass-im
//   (selects Dinv + esnid + output), --n-t0, --ninter, --gpu.
// Inputs : data_<free|ensbase>_Kdense_L<L>/K.<k>.h5   (mass-free; datasets <proj>_t<t>)
//          data_<ESNID>/prop_deter_L<L>/Dinv.<k>.h5    (ESNID = <...>_vmRe<re>vmIm<im>)
// Output : data_<ESNID>/corr_deter_exact_L<L>/corr.<k>.h5    (jj-style keys; nhits=1; atomic write).
// NOTE: vector tp/sp; ylm tower + axial are the next increment.
// -----------------------------------------------------------------------------

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <memory>
#include <getopt.h>
#include <highfive/H5File.hpp>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

using Idx = std::int32_t;
using Complex = std::complex<double>;
using CuC = cuDoubleComplex;

#define CUDA_CHECK(x) do{ cudaError_t e=(x); if(e!=cudaSuccess){ \
  fprintf(stderr,"CUDA %s:%d %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); std::exit(1);} }while(0)

namespace Comp{
#ifndef N_REFINE_CLI
#define N_REFINE_CLI 1
#endif
  constexpr int N_REFINE = N_REFINE_CLI;
  constexpr int NS = 2;
  constexpr int Nt = 128;
  constexpr Idx N_SITES = 10*N_REFINE*N_REFINE + 2;
  constexpr Idx Nx = NS*N_SITES;
  constexpr Idx N = Nx*Nt;
}


// time-shift a flat index by dt (time is the OUTER block of size Nx)
static inline Idx tshift(Idx i, int dt){
  const Idx Nx = Comp::Nx;
  const int Nt = Comp::Nt;
  const int t = (int)(i/Nx);
  const Idx x = i%Nx;
  const int tt = ((t+dt)%Nt + Nt)%Nt;
  return (Idx)tt*Nx + x;
}

// ANTIPERIODIC time BC: shifting a timeslice t by dt picks up (-1) for each wrap past the temporal
// boundary.  Returns the wrap sign (-1)^{#wraps} for moving time t by dt.
static inline int twrap_sign(int t, int dt){
  const int Nt = Comp::Nt;
  int s = t + dt;
  int w = 0;
  while(s < 0){ s += Nt; w++; }
  while(s >= Nt){ s -= Nt; w++; }
  return (w % 2 == 0) ? 1 : -1;
}


// ---- read a length-N^2 (or any) complex matrix from key/{real,imag} ----
static void load_mat(HighFive::File& f, const std::string& key, std::vector<Complex>& M){
  std::vector<double> re, im;
  f.getDataSet(key + "/real").read(re);
  f.getDataSet(key + "/imag").read(im);
  M.resize(re.size());
  for(size_t i=0; i<re.size(); i++) M[i] = Complex(re[i], im[i]);
}

// ---- write a length-Nt correlator: g = (1/4pi) C[t]; im negated for the conj channel ----
static void write_corr(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C, bool conj){
  const double inv4pi = 1.0/(4.0*M_PI);
  std::vector<double> re(Comp::Nt), im(Comp::Nt);
  for(int t=0; t<Comp::Nt; t++){
    const Complex g = inv4pi*C[t];
    re[t] = g.real();
    im[t] = conj ? -g.imag() : g.imag();
  }
  h5.createDataSet(key + "/real", re);
  h5.createDataSet(key + "/imag", im);
}

// ---- write a length-Nt complex vector RAW (no 1/4pi) ----
static void write_vec(HighFive::File& h5, const std::string& key, const std::vector<Complex>& C){
  std::vector<double> re(Comp::Nt), im(Comp::Nt);
  for(int t=0; t<Comp::Nt; t++){ re[t] = C[t].real(); im[t] = C[t].imag(); }
  h5.createDataSet(key + "/real", re);
  h5.createDataSet(key + "/imag", im);
}

// ---- A = K . P  (both row-major) via cuBLAS Zgemm. ----
// Trick: passing (P,K) to a column-major gemm computes (P^T)(K^T) = (K P)^T = the row-major K P.
static void matmul_A(cublasHandle_t cub, CuC* d_K, CuC* d_P, CuC* d_A,
                     const std::vector<Complex>& K, const std::vector<Complex>& P, std::vector<Complex>& A){
  const int n = (int)Comp::N;
  const CuC one  = make_cuDoubleComplex(1.0, 0.0);
  const CuC zero = make_cuDoubleComplex(0.0, 0.0);
  CUDA_CHECK(cudaMemcpy(d_K, reinterpret_cast<const CuC*>(K.data()), (size_t)n*n*sizeof(CuC), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_P, reinterpret_cast<const CuC*>(P.data()), (size_t)n*n*sizeof(CuC), cudaMemcpyHostToDevice));
  cublasZgemm(cub, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_P, n, d_K, n, &zero, d_A, n);
  A.resize((size_t)n*n);
  CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(A.data()), d_A, (size_t)n*n*sizeof(CuC), cudaMemcpyDeviceToHost));
}

// ---- tr(A) ----
static Complex trace(const std::vector<Complex>& A){
  const Idx N = Comp::N;
  Complex s(0,0);
  for(Idx i=0; i<N; i++) s += A[(size_t)i*N + i];
  return s;
}

// ---- tr(A B) ----
static Complex tr_AB(const std::vector<Complex>& A, const std::vector<Complex>& B){
  const Idx N = Comp::N;
  Complex s(0,0);
  for(Idx i=0; i<N; i++)
    for(Idx j=0; j<N; j++) s += A[(size_t)i*N + j] * B[(size_t)j*N + i];
  return s;
}

// ---- conn(dt) = tr(A . S_dt A S_-dt)  (free-field time-shift, ANTIPERIODIC) ----
// A(dt)[q,p] = (S_dt A S_-dt)[q,p] = eps_q eps_p A[tshift(q,-dt), tshift(p,-dt)], where eps = the
// antiperiodic wrap sign (twrap_sign) for shifting that index's time by -dt.  For the free
// (time-translation-invariant) field this reconstructs A(dt) from A0 EXACTLY.  Precompute the per-
// timeslice sign once per dt.
static Complex conn_shift(const std::vector<Complex>& A, int dt){
  const Idx N = Comp::N;
  const Idx Nx = Comp::Nx;
  const int Nt = Comp::Nt;
  std::vector<int> sgn(Nt);
  for(int t=0; t<Nt; t++) sgn[t] = twrap_sign(t, -dt);
  Complex s(0,0);
  for(Idx p=0; p<N; p++){
    const int sp = sgn[(int)(p/Nx)];
    for(Idx q=0; q<N; q++){
      const int sq = sgn[(int)(q/Nx)];
      // s += A[p,q] * A(dt)[q,p] ; eps_q eps_p = sq*sp
      s += A[(size_t)p*N + q] * (double)(sp*sq) * A[(size_t)tshift(q,-dt)*N + tshift(p,-dt)];
    }
  }
  return s;
}


// ---- CLI ----
struct Args {
  double mass_re = 0.0, mass_im = 0.0;
  std::string ens_dir;
  int n_t0 = 2, ninter = 10, gpu = 0;
};

void PrintHelp(){
  printf("jj_contract_deter: --mass-re --mass-im --ens-dir(omit=free) --n-t0 --ninter --gpu\n");
}

void ParseArgs(int argc, char** argv, Args& a){
  static struct option lo[] = {
    {"mass-re", required_argument, 0, 'r'},
    {"mass-im", required_argument, 0, 'i'},
    {"ens-dir", required_argument, 0, 'e'},
    {"n-t0",    required_argument, 0, 'T'},
    {"ninter",  required_argument, 0, 'I'},
    {"gpu",     required_argument, 0, 'G'},
    {"help",    no_argument,       0, 'h'},
    {0,0,0,0}
  };
  int opt, idx;
  while((opt = getopt_long(argc, argv, "r:i:e:T:I:G:h", lo, &idx)) != -1){
    switch(opt){
      case 'r': a.mass_re = std::stod(optarg); break;
      case 'i': a.mass_im = std::stod(optarg); break;
      case 'e': a.ens_dir = optarg;            break;
      case 'T': a.n_t0    = std::stoi(optarg); break;
      case 'I': a.ninter  = std::stoi(optarg); break;
      case 'G': a.gpu     = std::stoi(optarg); break;
      case 'h': default: PrintHelp(); std::exit(0);
    }
  }
}


int main(int argc, char* argv[]){
  std::cout << std::scientific << std::setprecision(15);
  Args a;
  ParseArgs(argc, argv, a);
  (void)a.gpu;
  CUDA_CHECK(cudaSetDevice(0));

  constexpr Idx N  = Comp::N;
  constexpr int Nt = Comp::Nt;
  const bool free_field = a.ens_dir.empty();
  const bool parity = (a.mass_im != 0.0);
  std::cout << "# contract: mass=(" << a.mass_re << "," << a.mass_im << ")  N=" << N
            << (free_field ? "  [free: time-shift]" : "  [interacting: per-t K]") << std::endl;

  // ---- output / input dirs (mirror the propagator + kbuild naming) ----
  std::string ens_base = a.ens_dir;
  if(!ens_base.empty() && ens_base.back()=='/') ens_base.pop_back();
  { const auto s = ens_base.find_last_of('/'); if(s!=std::string::npos) ens_base = ens_base.substr(s+1); }
  const std::string tag   = (free_field ? std::string("free") : ens_base);
  const std::string esnid = tag + "_vmRe" + std::to_string(a.mass_re) + "vmIm" + std::to_string(a.mass_im);
  const std::string kdir    = "data_" + tag   + "_Kdense_L"   + std::to_string(Comp::N_REFINE) + "/";  // mass-free
  const std::string propdir = "data_" + esnid + "/prop_deter_L" + std::to_string(Comp::N_REFINE) + "/";
  const std::string outdir  = "data_" + esnid + "/corr_deter_exact_L" + std::to_string(Comp::N_REFINE) + "/";
  std::filesystem::create_directories(outdir);

  const int n_t0 = a.n_t0;
  std::vector<int> t0s(n_t0);
  for(int b=0; b<n_t0; b++) t0s[b] = b*(Nt/n_t0);

  // ---- cuBLAS + device buffers for A = K . P ----
  cublasHandle_t cub;
  cublasCreate(&cub);
  CuC *d_K=nullptr, *d_P=nullptr, *d_A=nullptr;
  CUDA_CHECK(cudaMalloc(&d_K, (size_t)N*N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_P, (size_t)N*N*sizeof(CuC)));
  CUDA_CHECK(cudaMalloc(&d_A, (size_t)N*N*sizeof(CuC)));

  // ---- config loop ----
  const int k_lo   = 0;
  const int k_step = free_field ? 1 : a.ninter;
  const int k_hi   = free_field ? 1 : 1000000;
  for(int k=k_lo; k<k_hi; k+=k_step){
    const std::string kfile = kdir + "K." + std::to_string(k) + ".h5";
    const std::string pfile = propdir + "Dinv." + std::to_string(k) + ".h5";
    if(!std::filesystem::exists(kfile) || !std::filesystem::exists(pfile)){
      std::cout << "# k=" << k << " missing K/Dinv -- run stages 1a/1b first\n";
      if(free_field) break; else continue;
    }
    const std::string h5path = outdir + "corr." + std::to_string(k) + ".h5";
    if(std::filesystem::exists(h5path)){
      bool c=false;
      try { HighFive::File f(h5path, HighFive::File::ReadOnly); c = f.exist("complete"); } catch(...) {}
      if(c){ std::cout << "# skip k=" << k << " (complete)\n"; if(free_field) break; else continue; }
    }

    std::vector<Complex> P;
    { HighFive::File f(pfile, HighFive::File::ReadOnly); load_mat(f, "Dm_inv", P); }
    HighFive::File kf(kfile, HighFive::File::ReadOnly);   // open K once; datasets <proj>_t<t>
    std::cout << "# k=" << k << "  loaded P + opened K\n";

    // ATOMIC output: write to .tmp, rename after the 'complete' sentinel + close.
    const std::string h5tmp = h5path + ".tmp";
    auto h5p = std::make_unique<HighFive::File>(h5tmp,
                 HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate);
    HighFive::File& h5 = *h5p;
    h5.createDataSet("t0s",   t0s);
    h5.createDataSet("n_t0",  std::vector<int>{n_t0});
    h5.createDataSet("nhits", std::vector<int>{1});
    h5.createDataSet("ls",    std::vector<int>{0,1,2});

    for(const std::string proj : {std::string("tp"), std::string("sp")}){
      std::vector<std::vector<Complex>> Cpp(n_t0, std::vector<Complex>(Nt, Complex(0,0)));
      std::vector<Complex> discvec(Nt, Complex(0,0));

      if(free_field){
        // build A0 = K^P(0) P once; A(t) = S_t A0 S_-t => conn(dt) via the time-shift; disc = tr(A0).
        std::vector<Complex> K0, A0;
        load_mat(kf, proj + "_t0", K0);
        matmul_A(cub, d_K, d_P, d_A, K0, P, A0);
        const Complex disc = trace(A0);
        for(int t=0; t<Nt; t++) discvec[t] = disc;
        std::vector<Complex> cdt(Nt);
        for(int dt=0; dt<Nt; dt++) cdt[dt] = conn_shift(A0, dt);
        for(int b=0; b<n_t0; b++)
          for(int t=0; t<Nt; t++){ const int dt = ((t-t0s[b])%Nt + Nt)%Nt; Cpp[b][dt] = cdt[dt]; }
      } else {
        // interacting: no shift -> hold A(t0_b), stream A(t); conn = tr(A(t0) A(t)).
        std::vector<std::vector<Complex>> Ab(n_t0);
        for(int b=0; b<n_t0; b++){
          std::vector<Complex> Kb;
          load_mat(kf, proj + "_t" + std::to_string(t0s[b]), Kb);
          matmul_A(cub, d_K, d_P, d_A, Kb, P, Ab[b]);
        }
        for(int t=0; t<Nt; t++){
          std::vector<Complex> Kt, At;
          load_mat(kf, proj + "_t" + std::to_string(t), Kt);
          matmul_A(cub, d_K, d_P, d_A, Kt, P, At);
          discvec[t] = trace(At);
          for(int b=0; b<n_t0; b++){ const int dt = ((t-t0s[b])%Nt + Nt)%Nt; Cpp[b][dt] = tr_AB(Ab[b], At); }
        }
      }

      for(int b=0; b<n_t0; b++){
        const std::string kp = "h0/t0_" + std::to_string(b) + "/";
        write_corr(h5, kp + proj + "/Vpp", Cpp[b], false);
        if(!parity) write_corr(h5, kp + proj + "/Vmm", Cpp[b], true);   // massless/m_F: Vmm = conj(Vpp)
      }
      write_vec(h5, "h0/disc/" + proj + "/J", discvec);                 // RAW (matches jj disc/<proj>/J)
      std::cout << "#   " << proj << ": disc(0)=(" << discvec[0].real() << "," << discvec[0].imag()
                << ")  conn(dt=4)=" << Cpp[0][4].real() << "\n";
    }

    h5.createDataSet("complete", std::vector<int>{1});
    h5p.reset();
    std::filesystem::rename(h5tmp, h5path);
    std::cout << "# wrote " << h5path << "\n";
    if(free_field) break;
  }

  cublasDestroy(cub);
  CUDA_CHECK(cudaFree(d_K));
  CUDA_CHECK(cudaFree(d_P));
  CUDA_CHECK(cudaFree(d_A));
  return 0;
}
