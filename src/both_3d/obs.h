#pragma once

template<typename Gauge>
Complex get_polyakov( const Gauge& U, const Idx ix ){
  double res=0.0;
  for(int s=0; s<Comp::Nt; s++){
    res += U.tp( s, ix );
  }
  return std::exp( Complex(0.0, 1.0)*res );
}

// template<class T>
// void write( const std::string& path, const std::vector<T>& vector ){
//   std::ofstream ofs(path);
//   ofs << std::scientific << std::setprecision(15);
//   ofs << vector;
// }

// template<class T>
// void read( const std::string& path, std::vector<T>& vector ){
//   std::ifstream file(path);
//   std::string str;
//   while (std::getline(file, str)){
//     std::istringstream iss(str);
//     T v; iss >> v;
//     vector.push_back( v );
//   }
// }







// struct Scalar {
//   Double v;

//   Scalar()
//     : v(0.0)
//   {}

//   Scalar( const Double v_ )
//     : v(v_)
//   {}

//   Scalar( const Scalar& other )
//     : v(other.v)
//   {}

//   void clear(){ v = 0.0; }

//   Scalar& operator+=(const Scalar& rhs){
//     v += rhs.v;
//     return *this;
//   }

//   Scalar& operator+=(const Double& rhs){
//     v += rhs;
//     return *this;
//   }

//   Scalar& operator/=(const Double& rhs){
//     v /= rhs;
//     return *this;
//   }

//   std::string print() const {
//     std::stringstream ss;
//     ss << std::scientific << std::setprecision(15);
//     ss << v;
//     return ss.str();
//   }

//   void print(std::FILE* stream) const {
//     fprintf( stream, "%0.15Le\t", v );
//   }



// };



// struct Corr {
//   std::vector<Double> v;

//   // Corr()
//   //   : v(Lx*Ly, 0.0)
//   // {}

//   Corr( const std::vector<Double> v_ )
//     : v(v_)
//   {}

//   Corr( const Corr& other )
//     : v(other.v)
//   {}

//   // Double& operator()(const Idx x, const Idx y) { return v[idx(x,y)]; }
//   // Double operator()(const Idx x, const Idx y) const { return v[idx(x,y)]; }

//   void clear(){ for(Idx i=0; i<v.size(); i++) v[i] = 0.0; }

//   Corr& operator+=(const Corr& rhs)
//   {
//     // assert( rhs.v.size()==Lx*Ly );
//     for(Idx i=0; i<v.size(); i++) v[i] += rhs.v[i];
//     return *this;
//   }

//   Corr& operator*=(const double& a)
//   {
//     for(Idx i=0; i<v.size(); i++) v[i] *= a;
//     return *this;
//   }

//   Corr& operator+=(const double& a)
//   {
//     for(Idx i=0; i<v.size(); i++) v[i] += a;
//     return *this;
//   }


//   Corr& operator+=(const std::vector<Double>& rhs)
//   {
//     // assert( rhs.size()==Lx*Ly );
//     for(Idx i=0; i<v.size(); i++) v[i] += rhs[i];
//     return *this;
//   }

//   Corr& operator/=(const Double& rhs)
//   {
//     for(Idx i=0; i<v.size(); i++) v[i] /= rhs;
//     return *this;
//   }


//   // void print() const {
//   //   for(int y=0; y<Ly; y++){
//   //     for(int x=0; x<Lx; x++) {
//   //       printf( "%0.15e\t", (*this)(x, y) );
//   //     }
//   //     printf("\n");
//   //   }
//   // }

//   // void print(std::FILE* stream) const {
//   //   for(Idx y=0; y<Ly; y++){
//   //     for(Idx x=0; x<Lx; x++) {
//   //       if( !is_site(x,y) ) continue;
//   //       fprintf( stream, "%0.15Le\t", (*this)(x, y) );
//   //     }
//   //     fprintf( stream, "\n");
//   //   }
//   // }

//   // std::string print() const {
//   //   std::stringstream ss;
//   //   ss << std::scientific << std::setprecision(15);
//   //   for(Idx y=0; y<Ly; y++){
//   //     for(Idx x=0; x<Lx; x++) {
//   //       if( !is_site(x,y) ) continue;
//   //       ss << (*this)(x, y) << " ";
//   //     }
//   //     ss << std::endl;
//   //   }
//   //   return ss.str();
//   // }

// };


// template<typename T> // T needs to have: .clear, +=, /= defined
// struct Obs {
//   std::string description;
//   int N;
//   std::function<T(const Spin&)> f;

//   T sum;
//   int counter;

//   // Obs() = delete;

//   Obs
//   (
//    const std::string& description_,
//    const int N_,
//    const std::function<T(const Spin&)>& f_
//    )
//     : description(description_)
//     , N(N_)
//     , f(f_)
//     , sum()
//     , counter(0)
//   {}

//   void clear(){
//     sum.clear();
//     counter = 0;
//   }

//   void meas( const Spin& s ) {
//     sum += f(s);
//     counter++;
//   }

//   // T mean() const {
//   //   T tmp(sum);
//   //   tmp /= counter;
//   //   return mean;
//   // }

//   void write_and_clear( const std::string& dir, const int label ){
//     const std::string filename = dir + description + "_" + std::to_string(label) + ".dat";

//     // std::ofstream of( filename, std::ios::out | std::ios::trunc );
//     // if(!of) assert(false);
//     // of << std::scientific << std::setprecision(15);
//     // of << sum.print();
//     // of.close();

//     FILE *stream = fopen(filename.c_str(), "w");
//     if (stream == NULL) assert(false);
//     std::ofstream of( filename );
//     sum /= counter;
//     sum.print( stream );
//     fclose( stream );

//     clear();
//   }

// };
