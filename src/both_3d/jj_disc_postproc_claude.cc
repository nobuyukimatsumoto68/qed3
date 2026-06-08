// jj_disc_postproc_claude.cc
// -----------------------------------------------------------------------------
// Post-process the DISCONNECTED single-current traces J(t) stored in corr.<k>.h5
// (by jj_corr_claude.cu) into the projected disconnected two-point
//
//     C_disc^P(dt) = (1/4pi) * (1/Nt) sum_{t0} Jbar^P(t0) Jbar^P((t0+dt) mod Nt),
//
// where Jbar^P(t) = (1/nhits) sum_h J_h^P(t) is the hit-averaged single-current
// trace for projection P in {tp, sp, ylm l}.  PDF v2-10 Eqs. (3.65) [disc diagram
// = product of single-current traces] and (3.39) [physical correlator
// = -C_conn + C_disc; the sum is done DOWNSTREAM in the notebook].
//
// This program reads ONLY the h{h}/disc/... keys; it never touches the connected
// h{h}/t0_b/... keys and never modifies corr.<k>.h5.  Output is a sibling file
// disc_proced.<k>.h5 in the SAME corr_nt0<N>_nhits<H>/ directory, with leaves
//   tp/{Vpp,Vmm}   sp/{Vpp,Vmm}   ylm/{Vpp,Vmm}/l{l}    (length Nt, dt-indexed),
// mirroring the connected leaf names + write_corr/write_corr_conj re/im convention
// (one 1/4pi; im negated for the conj channel) so the notebook reuses its loader.
//
// Channel pairing (PDF Eqs. 3.64-3.67):
//   Vpp = Jbar * Jbar          (from disc/*/J)
//   Vmm = conj(Vpp)            (non-parity: massless / m_F, mirrors Vmm = conj(Vpp))
//   Vmm = Jtilbar * Jtilbar    (parity m_P: from disc/*/Jtil)
//
// Estimator (decision 1): hit-averaged product (O(1/nhits) self-contraction bias,
// fine at large nhits).  Unbiased cross-hit variant is a one-line switch (see below).
//
// Build (no CUDA): g++ -std=c++17 jj_disc_postproc_claude.cc -lhdf5 -o jj_disc_postproc_claude.o
//   (HighFive is header-only; same include path as the .cu files.)
// Usage: ./jj_disc_postproc_claude.o <corr_dir>   e.g. data_free_.../corr_nt02_nhits4/
// -----------------------------------------------------------------------------

#include <highfive/H5File.hpp>
#include <complex>
#include <vector>
#include <string>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <regex>

using Complex = std::complex<double>;
namespace fs = std::filesystem;

static const double INV4PI = 1.0 / (4.0 * std::acos(-1.0));

// ---- read a length-Nt complex vector from key/real + key/imag ----------------
static std::vector<Complex> read_vec(HighFive::File& f, const std::string& key){
  std::vector<double> re, im;
  f.getDataSet(key + "/real").read(re);
  f.getDataSet(key + "/imag").read(im);
  std::vector<Complex> C(re.size());
  for(size_t t = 0; t < re.size(); t++) C[t] = Complex(re[t], im[t]);
  return C;
}

// ---- hit-average a disc leaf "disc/<sub>" over all h groups -> Jbar(t) --------
// sub e.g. "tp/J", "sp/J", "ylm/l1/J", "tp/Jtil".  Jbar = (1/nhits) sum_h J_h.
static std::vector<Complex> hit_average(HighFive::File& f, int nhits, const std::string& sub){
  std::vector<Complex> Jbar;
  for(int h = 0; h < nhits; h++){
    auto Jh = read_vec(f, "h" + std::to_string(h) + "/disc/" + sub);
    if(Jbar.empty()) Jbar.assign(Jh.size(), Complex(0, 0));
    for(size_t t = 0; t < Jh.size(); t++) Jbar[t] += Jh[t];
  }
  for(auto& z : Jbar) z /= double(nhits);
  return Jbar;
  // Unbiased cross-hit alternative (drops the h=h' diagonal; needs nhits>=2):
  //   form S(t)=sum_h J_h(t) and Q(t1,t2)=sum_h J_h(t1)J_h(t2), then in two_point()
  //   use [S(t0)S(t0+dt) - Q(t0,t0+dt)] / (nhits*(nhits-1)).  Not enabled here.
}

// ---- circular translation-averaged two-point: ------------------------------
//   C(dt) = (1/Nt) sum_{t0} A(t0) B((t0+dt) mod Nt)
static std::vector<Complex> two_point(const std::vector<Complex>& A,
                                      const std::vector<Complex>& B){
  const int Nt = (int)A.size();
  std::vector<Complex> C(Nt, Complex(0, 0));
  for(int dt = 0; dt < Nt; dt++){
    Complex s(0, 0);
    for(int t0 = 0; t0 < Nt; t0++) s += A[t0] * B[(t0 + dt) % Nt];
    C[dt] = s / double(Nt);
  }
  return C;
}

// ---- write a length-Nt correlator leaf, mirroring jj_corr write_corr ----------
// g = (1/4pi) * C[t];  re = Re(g);  im = conj ? -Im(g) : Im(g).
static void write_corr(HighFive::File& o, const std::string& key,
                       const std::vector<Complex>& C, bool conj){
  const int Nt = (int)C.size();
  std::vector<double> re(Nt), im(Nt);
  for(int t = 0; t < Nt; t++){
    const Complex g = INV4PI * C[t];
    re[t] = g.real();
    im[t] = conj ? -g.imag() : g.imag();
  }
  o.createDataSet(key + "/real", re);
  o.createDataSet(key + "/imag", im);
}

// ---- process one corr.<k>.h5 -> disc_proced.<k>.h5 ---------------------------
static bool process_file(const std::string& in_path, const std::string& out_path){
  // input must be complete (the measurement sentinel) -- partial files lack disc traces.
  {
    HighFive::File f(in_path, HighFive::File::ReadOnly);
    if(!f.exist("complete")){
      std::cout << "#   skip (input INCOMPLETE): " << in_path << std::endl;
      return false;
    }
  }
  // resume: skip if the output already carries its own "complete" sentinel.
  if(fs::exists(out_path)){
    try {
      HighFive::File o(out_path, HighFive::File::ReadOnly);
      if(o.exist("complete")){ std::cout << "#   skip (done): " << out_path << std::endl; return false; }
    } catch(...) {}
  }

  HighFive::File f(in_path, HighFive::File::ReadOnly);

  std::vector<int> nh_v;  f.getDataSet("nhits").read(nh_v);
  const int nhits = nh_v.empty() ? 1 : nh_v[0];
  std::vector<int> ls;    f.getDataSet("ls").read(ls);   // l values for the ylm tower (e.g. 0,1,2)

  // parity: dagger-leg tilde traces are present only in the parity (imaginary mass) run.
  const bool parity = f.exist("h0/disc/tp/Jtil");

  // hit-averaged single-current traces, per projection.
  const auto Jtp = hit_average(f, nhits, "tp/J");
  const auto Jsp = hit_average(f, nhits, "sp/J");
  std::vector<std::vector<Complex>> Jyl(ls.size());
  for(size_t i = 0; i < ls.size(); i++)
    Jyl[i] = hit_average(f, nhits, "ylm/l" + std::to_string(ls[i]) + "/J");

  std::vector<Complex> JtpT, JspT;
  std::vector<std::vector<Complex>> JylT(ls.size());
  if(parity){
    JtpT = hit_average(f, nhits, "tp/Jtil");
    JspT = hit_average(f, nhits, "sp/Jtil");
    for(size_t i = 0; i < ls.size(); i++)
      JylT[i] = hit_average(f, nhits, "ylm/l" + std::to_string(ls[i]) + "/Jtil");
  }

  // disconnected two-points, both channels.
  const auto Ctp = two_point(Jtp, Jtp);
  const auto Csp = two_point(Jsp, Jsp);

  HighFive::File o(out_path,
      HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
  // carry-through scalars for the downstream loader.
  o.createDataSet("nhits",  std::vector<int>{nhits});
  o.createDataSet("ls",     ls);
  o.createDataSet("parity", std::vector<int>{parity ? 1 : 0});

  // tp
  write_corr(o, "tp/Vpp", Ctp, /*conj=*/false);
  if(parity) write_corr(o, "tp/Vmm", two_point(JtpT, JtpT), /*conj=*/false);
  else       write_corr(o, "tp/Vmm", Ctp,                   /*conj=*/true);   // Vmm = conj(Vpp)
  // sp
  write_corr(o, "sp/Vpp", Csp, /*conj=*/false);
  if(parity) write_corr(o, "sp/Vmm", two_point(JspT, JspT), /*conj=*/false);
  else       write_corr(o, "sp/Vmm", Csp,                   /*conj=*/true);
  // ylm tower
  for(size_t i = 0; i < ls.size(); i++){
    const std::string lk = "ylm/Vpp/l" + std::to_string(ls[i]);
    const std::string lm = "ylm/Vmm/l" + std::to_string(ls[i]);
    const auto Cl = two_point(Jyl[i], Jyl[i]);
    write_corr(o, lk, Cl, /*conj=*/false);
    if(parity) write_corr(o, lm, two_point(JylT[i], JylT[i]), /*conj=*/false);
    else       write_corr(o, lm, Cl,                          /*conj=*/true);
  }

  o.createDataSet("complete", std::vector<int>{1});   // sentinel, written LAST
  std::cout << "#   wrote " << out_path << (parity ? "  [parity]" : "") << std::endl;
  return true;
}

int main(int argc, char** argv){
  if(argc < 2){
    std::cerr << "Usage: " << argv[0] << " <corr_dir>   "
                 "(a data_<ESNID>/corr_nt0<N>_nhits<H>/ directory)\n";
    return 1;
  }
  const std::string corr_dir = argv[1];
  if(!fs::is_directory(corr_dir)){
    std::cerr << "ERROR: not a directory: " << corr_dir << std::endl;
    return 1;
  }

  // match corr.<k>.h5 (NOT disc_proced.*.h5); preserve <k> in the output name.
  const std::regex re_corr("^corr\\.([0-9]+)\\.h5$");
  int n_done = 0, n_seen = 0;
  for(const auto& ent : fs::directory_iterator(corr_dir)){
    if(!ent.is_regular_file()) continue;
    std::smatch m;
    const std::string fname = ent.path().filename().string();
    if(!std::regex_match(fname, m, re_corr)) continue;
    n_seen++;
    const std::string k = m[1].str();
    const std::string in_path  = ent.path().string();
    const std::string out_path = (fs::path(corr_dir) / ("disc_proced." + k + ".h5")).string();
    try {
      if(process_file(in_path, out_path)) n_done++;
    } catch(const std::exception& e){
      std::cerr << "#   ERROR on " << in_path << " : " << e.what() << std::endl;
    }
  }
  std::cout << "# done: " << n_done << " written / " << n_seen << " corr files in " << corr_dir << std::endl;
  return 0;
}
