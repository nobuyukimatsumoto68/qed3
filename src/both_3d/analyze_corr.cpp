#include <iostream>
#include <iomanip>
#include <filesystem>

#include <cmath>

#include <cassert>
#include <vector>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;


struct Corr {
  std::vector<double> data;

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

  double operator[](const int i) const { return data[i]; }
  double& operator[](const int i) { return data[i]; }

  Corr& operator+=(const Corr& rhs){
    if( data.size()==0 ) {
      data = rhs.data;
    }
    else{
      assert( data.size()==rhs.data.size() );
      for(int i=0; i<data.size(); i++) data[i] += rhs.data[i];
    }
    return *this;
  }

  Corr& operator-=(const Corr& rhs){
    if( data.size()==0 ) {
      data = rhs.data;
    }
    else{
      assert( data.size()==rhs.data.size() );
      for(int i=0; i<data.size(); i++) data[i] -= rhs.data[i];
    }
    return *this;
  }

  friend Corr operator-(Corr v, const Corr& w) { v -= w; return v; }

  Corr& operator*=(const Corr& rhs){
    assert( data.size()==rhs.data.size() );
    for(int i=0; i<data.size(); i++) data[i] *= rhs[i];
    return *this;
  }

  friend Corr operator*(Corr v, const Corr& w) { v *= w; return v; }

  Corr& operator*(const double& rhs){
    for(int i=0; i<data.size(); i++) data[i] *= rhs;
    return *this;
  }

  Corr& operator*=(const double& rhs){
    for(int i=0; i<data.size(); i++) data[i] *= rhs;
    return *this;
  }

  Corr& operator/=(const double& rhs){
    for(int i=0; i<data.size(); i++) data[i] /= rhs;
    return *this;
  }
};



template<typename T>
struct Jackknife {
  std::vector<T> binavgs;
  std::vector<T> jackavgs;
  T est;
  T var;

  int nbins;
  const int binsize;

  Jackknife( const int binsize )
    : binsize(binsize)
  {}

  void do_it(){
    nbins = binavgs.size();

    jackavgs.clear();
    jackavgs.resize( nbins );
    for(int i=0; i<nbins; i++){
      for(int j=0; j<nbins; j++){
        if(i==j) continue;
        jackavgs[i] += binavgs[j];
      }
      jackavgs[i] /= nbins-1.0;
    }

    est.clear();
    for(int i=0; i<nbins; i++){
      est += jackavgs[i];
    }
    est /= nbins;

    var.clear();
    for(int i=0; i<nbins; i++){
      T diff = jackavgs[i] - est;
      var += diff * diff;
    }
    var /= nbins;
    var *= nbins-1.0;
  }

};




int main(int argc, char **argv){

  std::string path = argv[1];
  std::string header = argv[2];
  int conf_min = atoi(argv[3]);
  int conf_max = atoi(argv[4]);
  int interval = atoi(argv[5]);
  int Nt = atoi(argv[6]);

  int binsize = atoi(argv[7]);
  // int Nt = atoi(argv[6]);

  std::cout << "# path = " << path << std::endl;
  std::cout << "# header = " << header << std::endl;
  std::cout << "# conf_min = " << conf_min << std::endl;
  std::cout << "# interval = " << interval << std::endl;
  std::cout << "# Nt = " << Nt << std::endl;
  std::cout << "# binsize = " << binsize << std::endl;

  // std::vector<fs::path> files;

  // std::vector<Corr> binavgs;
  Jackknife<Corr> jk(binsize);

  {
    int tmp=0;
    for(int conf=conf_min; conf<conf_max; conf+= interval){
      std::string file = path+"/"+header+std::to_string(conf)+".dat";
      // std::cout << "# file = " << file << std::endl;
      if( !fs::exists( file ) ) break;
      tmp = conf;
    }
    conf_max = tmp;
  }

  int nbins = (conf_max-conf_min)/interval/binsize;
  std::cout << "# nbins = " << nbins << std::endl;
  jk.binavgs.resize(nbins);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ib=0; ib<nbins; ib++){
    Corr binned;
    for(int conf = conf_min+ib*interval*binsize; conf<conf_min+(ib+1)*interval*binsize; conf+= interval){
      std::string file = path+"/"+header+std::to_string(conf)+".dat";
      if( !fs::exists( file ) ) assert(false);
      Corr corr(file);
      binned += corr;
    }
    binned /= binsize;
    jk.binavgs[ib] = binned;
  }

  jk.do_it();

  std::cout << "# test " << std::endl;

  for( int i=0; i<jk.est.size(); i++ ){
    std::cout << jk.est[i] << " " << std::sqrt(jk.var[i]) << std::endl;
  }

  // Iterate through the directory and add files to the vector
  // for (const auto& entry : fs::directory_iterator(path)) {
  //   if (fs::is_regular_file(entry)) {
  //     files.push_back(entry.path());
  //   }
  // }

  // // Sort the files alphabetically
  // std::sort(files.begin(), files.end());

  // // Print the sorted file names
  // for (const auto& file : files) {
  //   std::cout << file.filename().string() << std::endl;
  // }

  return 0;
}
