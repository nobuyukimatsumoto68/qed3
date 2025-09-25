
#include <iostream>
#include <iomanip>
#include <filesystem>

#include <cmath>

#include <cassert>
#include <vector>
#include <fstream>
#include <sstream>

#include <sstream>

namespace fs = std::filesystem;

#include <stdfloat>
using Real = std::float128_t;


struct Corr {
  std::vector<Real> data;

  Corr(){}

  Corr( const std::string& file ){
    data.clear();
    readfile( file );
  }

  std::size_t size() const { return data.size(); }

  void clear() {
    data.clear();
  }

  void readfile( const std::string& file ){
    std::string str;
    std::ifstream ifs(file);
    if(!ifs) assert(false);

    while (std::getline(ifs, str)){
      std::istringstream iss(str);
      double v1;
      iss >> v1;
      data.push_back( v1 );
    }
  }

  Real operator[](const std::size_t i) const { return data[i]; }
  Real& operator[](const std::size_t i) { return data[i]; }

  Corr& operator+=(const Corr& rhs){
    if( data.size()==0 ) {
      data = rhs.data;
    }
    else{
      // std::cout << "# debug. " << data.size() << " " << rhs.data.size() << std::endl;
      assert( data.size()==rhs.data.size() );
      for(std::size_t i=0; i<data.size(); i++) data[i] += rhs.data[i];
    }
    return *this;
  }

  Corr& operator-=(const Corr& rhs){
    if( data.size()==0 ) {
      data = rhs.data;
    }
    else{
      assert( data.size()==rhs.data.size() );
      for(std::size_t i=0; i<data.size(); i++) data[i] -= rhs.data[i];
    }
    return *this;
  }

  friend Corr operator-(Corr v, const Corr& w) { v -= w; return v; }

  Corr& operator*=(const Corr& rhs){
    assert( data.size()==rhs.data.size() );
    for(std::size_t i=0; i<data.size(); i++) data[i] *= rhs[i];
    return *this;
  }

  friend Corr operator*(Corr v, const Corr& w) { v *= w; return v; }

  Corr& operator*(const Real& rhs){
    for(std::size_t i=0; i<data.size(); i++) data[i] *= rhs;
    return *this;
  }

  Corr& operator*=(const Real& rhs){
    for(std::size_t i=0; i<data.size(); i++) data[i] *= rhs;
    return *this;
  }

  Corr& operator/=(const Real& rhs){
    for(std::size_t i=0; i<data.size(); i++) data[i] /= rhs;
    return *this;
  }
};



template<typename T>
struct Jackknife {
  std::vector<T> binavgs;
  std::vector<T> jackavgs;
  T est;
  T var;

  std::size_t nbins;
  const std::size_t binsize;

  Jackknife( const std::size_t binsize )
    : binsize(binsize)
  {}

  void do_it(){
    nbins = binavgs.size();

    jackavgs.clear();
    jackavgs.resize( nbins );
    for(std::size_t i=0; i<nbins; i++){
      for(std::size_t j=0; j<nbins; j++){
        if(i==j) continue;
        jackavgs[i] += binavgs[j];
      }
      jackavgs[i] /= nbins-1.0;
    }

    est.clear();
    for(std::size_t i=0; i<nbins; i++){
      est += jackavgs[i];
    }
    est /= nbins;

    var.clear();
    for(std::size_t i=0; i<nbins; i++){
      T diff = jackavgs[i] - est;
      var += diff * diff;
    }
    var /= nbins;
    var *= nbins-1.0;
  }


  void half_done(){
    nbins = jackavgs.size();

    est.clear();
    for(std::size_t i=0; i<nbins; i++){
      est += jackavgs[i];
    }
    est /= nbins;

    var.clear();
    for(std::size_t i=0; i<nbins; i++){
      T diff = jackavgs[i] - est;
      var += diff * diff;
    }
    var /= nbins;
    var *= nbins-1.0;
  }


};




int main(int argc, char **argv){

  std::string path = argv[1];
  std::cout << "# path = " << path << std::endl;
  std::string header = argv[2];
  std::cout << "# header = " << header << std::endl;
  std::size_t conf_min = atoi(argv[3]);
  std::cout << "# conf_min = " << conf_min << std::endl;
  std::size_t conf_max = atoi(argv[4]);
  std::size_t interval = atoi(argv[5]);
  std::cout << "# interval = " << interval << std::endl;
  std::size_t Nt = atoi(argv[6]);
  std::cout << "# Nt = " << Nt << std::endl;

  std::size_t binsize = atoi(argv[7]);
  std::cout << "# binsize = " << binsize << std::endl;
  std::size_t prefix_max = atoi(argv[8]);

  Real at = atof(argv[9]);
  std::string path_save = argv[10];


  Jackknife<Corr> jk(binsize);

  std::vector<std::size_t> nbins;

  {
    for(std::size_t prefix=1; prefix<=prefix_max; prefix++){
      std::size_t conf_max_loc = conf_min;
      // std::size_t tmp=0;
      for(std::size_t conf=conf_min; conf<conf_max; conf+= interval){
        std::string file = path+std::to_string(prefix)+"/"+header+std::to_string(conf)+".dat";
        // std::cout << "# file = " << file << std::endl;
        if( !fs::exists( file ) ) break;
        Corr corr(file);
        if( corr.size()!=Nt ) break;
        conf_max_loc = conf;
      }
      nbins.push_back( (conf_max_loc-conf_min)/interval/binsize ); // safety
    }
  }

  std::size_t nbins_tot = 0;
  for(const std::size_t elem : nbins) nbins_tot += elem;
  std::cout << "# nbins = " << nbins_tot << std::endl;
  jk.binavgs.resize(nbins_tot);


  for(std::size_t prefix=1; prefix<=prefix_max; prefix++){
    std::cout << "# prefix = " << prefix << std::endl;
    std::size_t accum=0;
    for(std::size_t pre=1; pre<prefix; pre++) accum += nbins[pre-1];
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(std::size_t ib=0; ib<nbins[prefix-1]; ib++){
      Corr binned;
      for(std::size_t conf = conf_min+ib*interval*binsize; conf<conf_min+(ib+1)*interval*binsize; conf+= interval){
        std::string file = path+std::to_string(prefix)+"/"+header+std::to_string(conf)+".dat";
        if( !fs::exists( file ) ) assert(false);
        Corr corr(file);
        if( corr.size()!=Nt ) assert(false);
        binned += corr;
      }
      binned /= binsize;
      jk.binavgs[accum + ib] = binned;
    }
  }

  std::cout << "# do it " << std::endl;

  jk.do_it();

  std::cout << "# test " << std::endl;
  {
    std::ofstream file(path_save+"_corr.dat", std::ios::trunc);
    for( std::size_t i=0; i<jk.est.size(); i++ ){
      file << (double)(jk.est[i]) << " " << (double)(std::sqrt(jk.var[i])) << std::endl;
    }
  }

  std::cout << "# test " << std::endl;
  Jackknife<Corr> jk_meff(binsize);
  jk_meff.jackavgs.resize(nbins_tot);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(std::size_t i=0; i<nbins_tot; i++) {
    jk_meff.jackavgs[i].data.resize(Nt);
    std::cout << "i = " << i << std::endl;
    jk_meff.jackavgs[i][0] = 0.0;
    for(int t=1; t<Nt; t++) {
      jk_meff.jackavgs[i][t] = - 1.0/at * std::log( jk.jackavgs[i][t]/jk.jackavgs[i][t-1] );
    }
  }
  std::cout << "# half_done " << std::endl;
  jk_meff.half_done();
  std::cout << "# half_done2 " << std::endl;
  {
    std::ofstream file(path_save+"_meff.dat", std::ios::trunc);
    for( std::size_t i=0; i<jk_meff.est.size(); i++ ){
      file << (double)(jk_meff.est[i]) << " " << (double)(std::sqrt(jk_meff.var[i])) << std::endl;
    }
  }


  return 0;
}
