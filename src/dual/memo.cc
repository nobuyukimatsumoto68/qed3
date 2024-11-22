// std::cout << "nn_oriented = " << std::endl;
// for(auto elem : nn_oriented[0]){
//   std::cout << elem << " ";
// }
// std::cout << std::endl;

// std::cout << "link_oriented = " << std::endl;
// for(auto elem : link_oriented[0]){
//   std::cout << elem << " ";
// }
// std::cout << std::endl;






  // std::map<std::pair<int, int>, double> omega;
  // std::map<std::pair<int, int>, double> alpha;
  std::unordered_map<Link, double> omega;
  std::unordered_map<Link, double> alpha;
  std::unordered_map<int, int> NM2EO;

  {
    std::ifstream file("omega.dat");

    std::string str;
    std::string file_contents;
    while (std::getline(file, str)){
      // std::cout << str << std::endl;
      std::istringstream iss(str);
      int i,j;
      double v;
      iss >> i;
      iss >> j;
      iss >> v;
      // omega.insert( std::make_pair( std::pair<int, int>(i,j), v ) );
      omega.insert( { Link(i,j,lattice.n_sites), v } );
      omega.insert( { Link(j,i,lattice.n_sites), -v } );
    }

    std::cout << "omega = " << std::endl;
    for(auto elem : omega){
      std::cout << elem.first[0] << " " << elem.first[1] << " " << elem.second << std::endl;
      // std::cout << omega[std::make_pair(elem.first.first, elem.first.second)] << " " << elem.second << std::endl;
    }

    std::cout << "size of omega = " << omega.size() << std::endl;
  }

  {
    std::ifstream file("alpha.dat");

    std::string str;
    std::string file_contents;
    while (std::getline(file, str)){
      // std::cout << str << std::endl;
      std::istringstream iss(str);
      int i,j;
      double v;
      iss >> i;
      iss >> j;
      iss >> v;
      // alpha.insert( std::make_pair( std::pair<int, int>(i,j), v ) );
      alpha.insert( {Link(i,j,lattice.n_sites), v} );
      // alpha.insert( std::make_pair( Link{j,i}, -v ) );
    }

    std::cout << "alpha = " << std::endl;
    for(auto elem : alpha){
      // std::cout << elem.first.first << " " << elem.first.second << " " << elem.second << std::endl;
      std::cout << elem.first[0] << " " << elem.first[1] << " " << elem.second << std::endl;
    }
  }


  // {
  //   std::cout << "size of omega = " << omega.size() << std::endl;
  //   for(auto elem : omega){
  //     // int ix = elem.first.first;
  //     // int iy = elem.first.second;
  //     int ix = elem.first[0];
  //     int iy = elem.first[1];

  //     // double alpha1 = alpha[std::make_pair(ix,iy)];
  //     // double alpha2 = alpha[std::make_pair(iy,ix)];
  //     // double omega12 = omega[std::make_pair(ix,iy)];
  //     double alpha1 = alpha[Link{ix,iy}];
  //     double alpha2 = alpha[Link{iy,ix}];
  //     double omega12 = omega[Link{ix,iy}];

  //     double diff = (alpha2 + M_PI + omega12) - alpha1;
  //     std::cout << "ix,iy = " << ix << ", " << iy << ", diff = " << Mod(diff) << std::endl;
  //   }
  // }

  {
    NM2EO.insert( {10, 0} );

    NM2EO.insert( { 3, 3} );
    NM2EO.insert( { 9, 5} );
    NM2EO.insert( { 1, 8} );
    NM2EO.insert( { 7, 9} );
    NM2EO.insert( { 8,11} );

    NM2EO.insert( { 6,10} );
    NM2EO.insert( { 5, 1} );
    NM2EO.insert( { 2, 2} );
    NM2EO.insert( {12, 4} );
    NM2EO.insert( { 4, 7} );

    NM2EO.insert( {11, 6} );
  }

// for(auto elem : NM2EO){
//   std::cout << elem.first << " to " << elem.second << std::endl;
// }

// std::unordered_map<std::pair<int, int>, double> omegaEO;
// std::unordered_map<std::pair<int, int>, double> alphaEO;
std::unordered_map<Link, double> omegaEO;
std::unordered_map<Link, double> alphaEO;
{
  for(auto elem : alpha){
    // int ix1 = elem.first.first;
    // int iy1 = elem.first.second;
    int ix1 = elem.first[0];
    int iy1 = elem.first[1];
    int ix2 = NM2EO[ix1];
    int iy2 = NM2EO[iy1];
    std::cout << ix1 << ", " << iy1
	      << " mapped to "
	      << ix2 << ", " << iy2
	      << std::endl;
    // alphaEO.insert( std::make_pair( std::pair<int, int>(ix2,iy2),
    // 				      alpha[std::pair<int, int>(ix1,iy1)] ) );
    alphaEO.insert( { Link(ix2,iy2,lattice.n_sites), alpha[Link(ix1,iy1,lattice.n_sites)] } );
  }

  for(auto elem : omega){
    // int ix1 = elem.first.first;
    // int iy1 = elem.first.second;
    int ix1 = elem.first[0];
    int iy1 = elem.first[1];
    int ix2 = NM2EO[ix1];
    int iy2 = NM2EO[iy1];
    // omegaEO.insert( std::make_pair( std::pair<int, int>(ix2,iy2),
    // 				      omega[std::pair<int, int>(ix1,iy1)] ) );
    omegaEO.insert( { Link(ix2,iy2,lattice.n_sites), omega[Link(ix1,iy1,lattice.n_sites)] } );
  }

  // std::cout << "omegaEO = " << std::endl;
  // for(auto elem : omegaEO){
  //   std::cout << elem.first[0] << " " << elem.first[1] << " " << elem.second << std::endl;
  // }
  // std::cout << "alphaEO = " << std::endl;
  // for(auto elem : alphaEO){
  //   std::cout << elem.first[0] << " " << elem.first[1] << " " << elem.second << std::endl;
  // }
}

  {
    int ix1 = 1;
    int iy1 = 5;

    std::cout << std::endl;
    std::cout << ix1 << ", " << iy1
	      << " mapped to "
	      << NM2EO[ix1] << ", " << NM2EO[iy1]
	      << std::endl;
    std::cout << "omegaNM = " << omega[Link(ix1,iy1,lattice.n_sites)] << std::endl;
    std::cout << "omegaEO = " << omegaEO[Link(NM2EO[ix1],NM2EO[iy1],lattice.n_sites)] << std::endl;
    std::cout << "alphaNM = " << alpha[Link(ix1,iy1,lattice.n_sites)] << std::endl;
    std::cout << "alphaEO = " << alphaEO[Link(NM2EO[ix1],NM2EO[iy1],lattice.n_sites)] << std::endl;
    std::cout << "alphaNM^-1 = " << alpha[Link(iy1,ix1,lattice.n_sites)] << std::endl;
    std::cout << "alphaEO^-1 = " << alphaEO[Link(NM2EO[iy1],NM2EO[ix1],lattice.n_sites)] << std::endl;
    std::cout << std::endl;
  }


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// EO 1, 5
// NM 2, 9
  {
    // for(auto elem : omega){
    {
      int ix = 10;
      int iy = 3;

      // auto elem = omega[std::make_pair(ix,iy)];

      // double alpha1 = alpha[std::make_pair(ix,iy)];
      // double alpha2 = alpha[std::make_pair(iy,ix)];
      // double omega12 = omega[std::make_pair(ix,iy)];
      double alpha1 = alpha[Link(ix,iy,lattice.n_sites)];
      double alpha2 = alpha[Link(iy,ix,lattice.n_sites)];
      double omega12 = omega[Link(ix,iy,lattice.n_sites)];

      double diff = (alpha2 + M_PI + omega12) - alpha1;
      
      std::cout << "diff NM = " << Mod(diff) << std::endl;
      std::cout << "alpha1 = " << alpha1 << std::endl;
      std::cout << "alpha2 = " << alpha2 << std::endl;
      std::cout << "omega12 = " << omega12 << std::endl;
    }

    {
      int ix = 1;
      int iy = 5;
     
      // auto elem = omegaEO[std::make_pair(ix,iy)];

      // double alpha1 = alphaEO[std::make_pair(ix,iy)];
      // double alpha2 = alphaEO[std::make_pair(iy,ix)];
      // double omega12 = omegaEO[std::make_pair(ix,iy)];
      double alpha1 = alphaEO[Link(ix,iy,lattice.n_sites)];
      double alpha2 = alphaEO[Link(iy,ix,lattice.n_sites)];
      double omega12 = omegaEO[Link(ix,iy,lattice.n_sites)];

      double diff = (alpha2 + M_PI + omega12) - alpha1;
      std::cout << "diff EO = " << Mod(diff) << std::endl;
      std::cout << "alpha1 = " << alpha1 << std::endl;
      std::cout << "alpha2 = " << alpha2 << std::endl;
      std::cout << "omega12 = " << omega12 << std::endl;
    }
  }
  


  {
    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<lattice.sites[ix].nn; jj++){
	const int iy = lattice.sites[ix].neighbors[jj];

	double alpha1 = alphaEO[Link(ix,iy,lattice.n_sites)];
	double alpha2 = alphaEO[Link(iy,ix,lattice.n_sites)];
	double omega12 = omegaEO[Link(ix,iy,lattice.n_sites)];

	double diff = (alpha2 + M_PI + omega12) - alpha1;
	std::cout << "ix,iy = " << ix << ", " << iy << ", diff = " << Mod(diff) << std::endl;
      }}
  }






  // return 0;

  // Eigen::IOFormat fmt(15, 0,
  // 		      "\t",
  // 		      "\n",
  // 		      "",
  // 		      "",
  // 		      "",
  // 		      "");

  // ----------------------------------

  // std::unordered_map< Link, double > omega;
  // map.insert( std::make_pair( Link(3,10) ) )
  // std::vector<double> rationals_NM = ;

  // ----------------------------------


  // std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);

  // // const int seed = 1;
  // // std::mt19937 gen(seed);
  // for(auto elem : D.alpha0){
  //   elem = 0.0; // dis(gen);
  // }
  // D.set_omega();

  // D.omega[0] = double_mod(D.omega[0] - M_PI/5.0);
  // D.omega[1] = double_mod(D.omega[0] - M_PI/5.0);

  // // D.alpha0[0] = M_PI/2.0;

  // // {
  // //   const int iN = 0;
  // //   const int iS = 6;
  // //   assert( std::abs(lattice.r[iN][2]-1.0)<1.0e-15 );
  // //   assert( std::abs(lattice.r[iS][2]+1.0)<1.0e-15 );

  // //   double alpha0_N = 0.0;
  // //   double alpha0_S = M_PI/2.0;

  // //   { // N
  // //     D.alpha[iN].resize(lattice.sites[iN].nn);
  // //     for(int jj=0; jj<lattice.sites[iN].nn; jj++){
  // // 	const int il = D.link_oriented[iN][jj];
  // // 	D.omega[il] = M_PI/5.0;

  // // 	D.alpha[iN][jj] = double_mod( alpha0_N + 2.0*M_PI/5.0 * jj );
  // //     }
  // //   } // end N
  // //   { // S
  // //     D.alpha[iS].resize(lattice.sites[iS].nn);
  // //     for(int jj=0; jj<lattice.sites[iS].nn; jj++){
  // // 	const int il = D.link_oriented[iS][jj];
  // // 	D.omega[il] = M_PI/5.0;

  // // 	D.alpha[iS][jj] = double_mod( alpha0_S + 2.0*M_PI/5.0 * jj );
  // //     }
  // //   } // end S
  // //   { // nn of N
  // //     for(int jj=0; jj<lattice.sites[iN].nn; jj++){
  // // 	const int il = D.link_oriented[iN][jj];
  // // 	const int ix = D.nn_oriented[iN][jj];
  // // 	D.alpha[ix].resize(lattice.sites[ix].nn);

  // // 	int kk0=0; // find kk0 that corresponds to il
  // // 	for(; kk0<lattice.sites[ix].nn; kk0++) if(D.link_oriented[ix][kk0]==il) break;
  // // 	assert(kk0!=lattice.sites[ix].nn);

  // // 	D.alpha[ix][kk0] = D.alpha[iN][jj] - D.omega[il];

  // // 	for(int dk=1; dk<lattice.sites[ix].nn; dk++){
  // // 	  const int kk = (kk0+dk)%lattice.sites[ix].nn;
  // // 	  D.alpha[ix][kk] = double_mod( D.alpha[ix][kk0] + 2.0*M_PI/5.0 * jj );
  // // 	}
  // //     }
  // //   } // end nn of N
  // //   { // nn of S
  // //     for(int jj=0; jj<lattice.sites[iS].nn; jj++){
  // // 	const int il = D.link_oriented[iS][jj];
  // // 	const int ix = D.nn_oriented[iS][jj];
  // // 	D.alpha[ix].resize(lattice.sites[ix].nn);

  // // 	int kk0=0; // find kk0 that corresponds to il
  // // 	for(; kk0<lattice.sites[ix].nn; kk0++) if(D.link_oriented[ix][kk0]==il) break;
  // // 	assert(kk0!=lattice.sites[ix].nn);
  // // 	D.alpha[ix][kk0] = D.alpha[iS][jj] - D.omega[il];

  // // 	for(int dk=1; dk<lattice.sites[ix].nn; dk++){
  // // 	  const int kk = (kk0+dk)%lattice.sites[ix].nn;
  // // 	  D.alpha[ix][kk] = double_mod( D.alpha[ix][kk0] + 2.0*M_PI/5.0 * jj );
  // // 	}
  // //     }
  // //   } // end nn of S
  // //   { // set omega
  // //     for(int il=0; il<lattice.n_links; il++){
  // // 	const QfeLink link = lattice.links[il];
  // // 	const int ix = std::min(link.sites[0], link.sites[1]);
  // // 	const int iy = std::max(link.sites[0], link.sites[1]);

  // // 	const auto itx_ell = std::find(D.link_oriented[ix].begin(),
  // // 				       D.link_oriented[ix].end(),
  // // 				       il);
  // // 	const auto ity_ell = std::find(D.link_oriented[iy].begin(),
  // // 				       D.link_oriented[iy].end(),
  // // 				       il);
  // // 	const int ixl = std::distance(D.link_oriented[ix].begin(), itx_ell);
  // // 	const int iyl = std::distance(D.link_oriented[iy].begin(), ity_ell);

  // // 	D.omega[il] = double_mod( D.alpha[ix][ixl] - D.alpha[iy][iyl] - M_PI );
  // //     }
  // //   }
  // // }

  // // {
  // //   for(auto elem : D.omega) {
  // //     std::cout << std::setw(25) << elem << " ";
  // //   }
  // //   std::cout << std::endl;
  // //   // std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // // }

  // {
  //   for(int ia=0; ia<lattice.n_faces; ia++){
  //     double delta_alpha_sum = 0.0;
  //     double omega_sum = 0.0;

  //     for(int i=0; i<3; i++){
  // 	int il = lattice.faces[ia].edges[i];
  // 	omega_sum += D.omega[il];

  // 	// ----------------------------

  // 	const QfeLink link = lattice.links[il];
  // 	const int ix = std::min(link.sites[0], link.sites[1]);
  // 	const int iy = std::max(link.sites[0], link.sites[1]);

  // 	const auto itx_ell = std::find(D.link_oriented[ix].begin(),
  // 				       D.link_oriented[ix].end(),
  // 				       il);
  // 	const auto ity_ell = std::find(D.link_oriented[iy].begin(),
  // 				       D.link_oriented[iy].end(),
  // 				       il);
  // 	const int ixl = std::distance(D.link_oriented[ix].begin(), itx_ell);
  // 	const int iyl = std::distance(D.link_oriented[iy].begin(), ity_ell);

  // 	// double alpha_x = D.alpha[ix][ixl];
  // 	// double alpha_y = D.alpha[iy][iyl] + M_PI;
  // 	double alpha_x = D.alpha(ix,ixl);
  // 	double alpha_y = D.alpha(iy,iyl) + M_PI;
  // 	double diff = alpha_x - (alpha_y+D.omega[il]);
  // 	std::cout << "diff = " << double_mod( diff ) << std::endl;

  // 	// ----------------------------
  //     }

  //     std::cout << "ia = " << ia
  // 		<< ", sum (directed) = " << double_mod( D.face_signs[ia]*omega_sum ) << std::endl;
      
  //     // ------------------------

  //   }
  // }



  // std::cout << "# x y z" << std::endl;
  // std::cout << std::scientific << std::setprecision(15);

  // {
  //   auto vec = lattice.r[ix];
  //   for(auto elem : vec) {
  //     std::cout << std::setw(25) << elem << " ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // }
  // {
  //   auto vec = lattice.r[iy];
  //   for(auto elem : vec) {
  //     std::cout << std::setw(25) << elem << " ";
  //   }
  //   std::cout << std::endl;
  //   std::cout << projection(vec).format(fmt).transpose() << std::endl;
  // }




	// int il = lattice.faces[ia].edges[i];


	// omega_sum += D.omega[il];

	// // ----------------------------

	// const QfeLink link = lattice.links[il];
	// const int ix = std::min(link.sites[0], link.sites[1]);
	// const int iy = std::max(link.sites[0], link.sites[1]);

	// const auto itx_ell = std::find(D.link_oriented[ix].begin(),
	// 			       D.link_oriented[ix].end(),
	// 			       il);
	// const auto ity_ell = std::find(D.link_oriented[iy].begin(),
	// 			       D.link_oriented[iy].end(),
	// 			       il);
	// const int ixl = std::distance(D.link_oriented[ix].begin(), itx_ell);
	// const int iyl = std::distance(D.link_oriented[iy].begin(), ity_ell);

	// // double alpha_x = D.alpha[ix][ixl];
	// // double alpha_y = D.alpha[iy][iyl] + M_PI;
	// double alpha_x = D.alpha(ix,ixl);
	// double alpha_y = D.alpha(iy,iyl) + M_PI;
	// double diff = alpha_x - (alpha_y+D.omega[il]);
	// std::cout << "diff = " << double_mod( diff ) << std::endl;
	// ----------------------------


  std::cout << SW.gR << std::endl;
  {
    std::cout << "vps = " << std::endl;
    int i=0;
    double sum = 0.0;
    for(auto elem : lattice.vps) {
      sum += elem;
      i++;
    }
    std::cout << "sum = " << sum << std::endl;
    std::cout << "check = " << 4.0*M_PI << std::endl;
  }
  std::cout << U.plaquette_angle( 0 ) << std::endl;
  std::cout << SW( U ) << std::endl;



      // if(is_compact) plaq.Measure( U.average_plaquette() );
      // else
      // double val = U.plaquette_angle(iface);
      // while(val>M_PI) val -= 2.0*M_PI;
      // while(val<-M_PI) val += 2.0*M_PI;
      // std::cout << val << std::endl;


  // std::cout << var.Mean() << " ";
  // std::cout << var.Error() * std::sqrt( var.AutocorrBack() )
  // 	    << std::endl;
  // }
  // << " "
  // << plaq.AutocorrFront() << " "
  // << plaq.AutocorrBack()

  
  // double denom = 0.0;
  // double numer = 0.0;

  // double tmp = 0.0;
  // // const int n=2;
  // // for(int n=0; n<lattice.n_faces; n++){
  // //   std::cout << std::exp(-beta) * boost::math::cyl_bessel_i(n, beta) << std::endl;
  // //   const double bessel_n = ine(n, beta);
  // //   const double bessel_np1 = ine(n+1, beta);
  // //   const double bessel_nm1 = ine(std::abs(n-1), beta);
  // //   const double d_bessel_n = 0.5*(bessel_np1 + bessel_nm1);
  // //   // std::cout << "bessel_n = " << bessel_n << std::endl;
  // //   // std::cout << "bessel_np1 (boost)" << std::exp(-beta) * boost::math::cyl_bessel_i(n+1, beta) << std::endl;
  // //   // std::cout << "bessel_np1 = " << bessel_np1 << std::endl;
  // //   // std::cout << "bessel_nm1 = " << bessel_nm1 << std::endl;
  // //   // std::cout << "d_bessel_n = " << d_bessel_n << std::endl;

  // //   denom += std::pow(bessel_n / bessel_0, lattice.n_faces);
  // //   numer += d_bessel_n / bessel_0 * std::pow(bessel_n / bessel_0, lattice.n_faces-1);

  // //   std::cout << "denom = " << denom << std::endl;
  // //   std::cout << "numer = " << numer << std::endl;
  // // }
  // double beta = 4.0;
  // const double bessel_0 = std::exp(-beta) * boost::math::cyl_bessel_i(0, beta);
  // for(int n=0; n<lattice.n_faces; n++){
  //   const double bessel_n = std::exp(-beta) * boost::math::cyl_bessel_i(n, beta);
  //   const double bessel_np1 = std::exp(-beta) * boost::math::cyl_bessel_i(n+1, beta);
  //   const double bessel_nm1 = std::exp(-beta) * boost::math::cyl_bessel_i(std::abs(n-1), beta);
  //   const double d_bessel_n = 0.5*(bessel_np1 + bessel_nm1);
  //   // const double bessel_n = ine(n, beta);
  //   // std::cout << "bessel = " << bessel_n << std::endl;
  //   // const double bessel_np1 = ine(n+1, beta);
  //   // std::cout << "bessel = " << bessel_np1 << std::endl;
  //   // const double bessel_nm1 = ine(std::abs(n-1), beta);
  //   // std::cout << "bessel = " << bessel_nm1 << std::endl;
  //   // const double d_bessel_n = 0.5*(bessel_np1 + bessel_nm1);

  //   denom += std::pow(bessel_n / bessel_0, lattice.n_faces);
  //   numer += d_bessel_n / bessel_0 * std::pow(bessel_n / bessel_0, lattice.n_faces-1);
  //   // std::cout << "denom = " << denom << std::endl;
  //   // std::cout << "numer = " << numer << std::endl;
  //   // if( std::abs(numer / denom - tmp ) < 1.0e-10  ) break;
  //   // tmp = numer / denom;
  // }
  // const double exact = numer / denom;
  // std::cout << "exact = " << lattice.n_faces * ( 1.0 - exact ) << std::endl;


  // const int iface = 1;
  // const double factor = 0.5 / U.info_p.vp(iface);
  // std::cout << "vp = " << U.info_p.vp(iface) << std::endl;
  // const double factor = 0.5;


// double ine( const int nu, const double z ){
//   double res = 0.0;
//   double tmp = 0.0;
//   if(std::abs(z)>20.0){
//     double pow = 1.0;
//     for(int n=0; n<100; n++){
//       const double coeff = std::tgamma( nu+n+0.5 ) / std::tgamma( n+1 ) / std::tgamma( nu-n+0.5 );
//       // std::cout << "coeff = " << coeff << std::endl;
//       res += std::pow( -1.0/(2.0*z), n ) * coeff;
//       // res += pow * coeff;
//       // pow *= -1.0/(2.0*z);
//       if( std::abs(res - tmp) < 1.0e-10  ) break;
//       tmp = res;
//     }
//     res /= std::sqrt(2.0*M_PI*z);
//   }
//   else{
//     res = std::exp(-z) * boost::math::cyl_bessel_i(nu, z);
//   }
//   return res;
// }




  { // first in coo
    const int N = 2*lattice.n_sites;
    std::vector<int> cols;
    std::vector<int> rows;

    // std::vector<int> passive_counter_cols_visited(N,0);
    // std::vector<int> passive_counter_rows_visited(N,0);
    // std::vector<int> active_counter_cols_visited;
    // std::vector<int> active_counter_rows_visited;

    // int counter=0;
    // rows.push_back(counter);
    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];

	rows.push_back( NS*ix ); cols.push_back( NS*iy );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*iy ] );
	// ++passive_counter_rows_visited[ Ns*ix ];
	// ++passive_counter_cols_visited[ Ns*iy ];

	rows.push_back( NS*ix ); cols.push_back( NS*iy+1 );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*iy+1 ] );
	// ++passive_counter_rows_visited[ Ns*ix ];
	// ++passive_counter_cols_visited[ Ns*iy+1 ];

	// res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
	rows.push_back( NS*ix ); cols.push_back( NS*ix );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*ix ] );
	// ++passive_counter_rows_visited[ Ns*ix ];
	// ++passive_counter_cols_visited[ Ns*ix ];

	rows.push_back( NS*ix ); cols.push_back( NS*ix+1 );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*ix+1 ] );
	// ++passive_counter_rows_visited[ Ns*ix ];
	// ++passive_counter_cols_visited[ Ns*ix+1 ];

	// res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
	rows.push_back( NS*ix+1 ); cols.push_back( NS*iy );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix+1 ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*iy ] );
	// ++passive_counter_rows_visited[ Ns*ix+1 ];
	// ++passive_counter_cols_visited[ Ns*iy ];

	rows.push_back( NS*ix+1 ); cols.push_back( NS*iy+1 );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix+1 ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*iy+1 ] );
	// ++passive_counter_rows_visited[ Ns*ix+1 ];
	// ++passive_counter_cols_visited[ Ns*iy+1 ];

	// res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
	rows.push_back( NS*ix+1 ); cols.push_back( NS*ix );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix+1 ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*ix ] );
	// ++passive_counter_rows_visited[ Ns*ix+1 ];
	// ++passive_counter_cols_visited[ Ns*ix ];

	rows.push_back( NS*ix+1 ); cols.push_back( NS*ix+1 );
	// active_counter_rows_visited.push_back( passive_counter_rows_visited[ Ns*ix+1 ] );
	// active_counter_cols_visited.push_back( passive_counter_cols_visited[ Ns*ix+1 ] );
	// ++passive_counter_rows_visited[ Ns*ix+1 ];
	// ++passive_counter_cols_visited[ Ns*ix+1 ];
      }
    }

    const int len = cols.size();

    std::vector<int> cols_csr(len,-1);
    std::vector<int> rows_csr(len,-1);

    std::vector<int> cols_csrT(len,-1);
    std::vector<int> rows_csrT(len,-1);

    int counter_csr=0;
    int counter_csrT=0;

    rows_csr( counter_csr );
    for(int i=0; i<N; i++){
      for(int k=0; k<len; k++){
	if( rows[k]==i ){
	  cols_csr[k] = counter_csr;
	  ++counter_csr;
	}
	if( cols[k]==i ){
	  cols_csrT[k] = counter_csrT;
	  ++counter_csrT;
	}
      }
      rows_csr( counter_csr );
      rows_csrT( counter_csrT );
    }
    
  }


  // ugliness of repetition
  void csr( std::vector<Complex>& v,
	    std::vector<int>& cols,
	    std::vector<int>& rows,
	    const U1onS2& U ) const {
    int counter=0;
    rows.push_back(counter);

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
	const MS tmp2 = lattice.vol[ix]/lattice.mean_vol * r*lattice.u[ix][jj]*sigma[0];

	// res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
	v.push_back( -tmp(0,0) ); cols.push_back( NS*iy ); counter++;
	v.push_back( -tmp(0,1) ); cols.push_back( NS*iy+1 ); counter++;

	// res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
	v.push_back( tmp2(0,0) ); cols.push_back( NS*ix ); counter++;
	v.push_back( tmp2(0,1) ); cols.push_back( NS*ix+1 ); counter++;
      }
      rows.push_back(counter);

      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
	const MS tmp2 = lattice.vol[ix]/lattice.mean_vol * r*lattice.u[ix][jj]*sigma[0];

	// res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
	v.push_back( -tmp(1,0) ); cols.push_back( NS*iy ); counter++;
	v.push_back( -tmp(1,1) ); cols.push_back( NS*iy+1 ); counter++;

	// res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
	v.push_back( tmp2(1,0) ); cols.push_back( NS*ix ); counter++;
	v.push_back( tmp2(1,1) ); cols.push_back( NS*ix+1 ); counter++;
      }
      rows.push_back(counter);
    }
  }


  void coo( std::vector<Complex>& v,
	    std::vector<int>& cols,
	    std::vector<int>& rows,
	    const U1onS2& U ) const {

    for(int ix=0; ix<lattice.n_sites; ix++){
      for(int jj=0; jj<3; jj++){
	int iy = lattice.nns[ix][jj];
	const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
	const MS tmp2 = lattice.vol[ix]/lattice.mean_vol * r*lattice.u[ix][jj]*sigma[0];

	// res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
	v.push_back( -tmp(0,0) ); rows.push_back( NS*ix ); cols.push_back( NS*iy );
	v.push_back( -tmp(0,1) ); rows.push_back( NS*ix ); cols.push_back( NS*iy+1 );

	// res[NS*ix] += tmp(0,0)*v[NS*ix] + tmp(0,1)*v[NS*ix+1];
	v.push_back( tmp2(0,0) ); rows.push_back( NS*ix ); cols.push_back( NS*ix );
	v.push_back( tmp2(0,1) ); rows.push_back( NS*ix ); cols.push_back( NS*ix+1 );

	// res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
	v.push_back( -tmp(1,0) ); rows.push_back( NS*ix+1 ); cols.push_back( NS*iy );
	v.push_back( -tmp(1,1) ); rows.push_back( NS*ix+1 ); cols.push_back( NS*iy+1 );

	// res[NS*ix+1] += tmp(1,0)*v[NS*ix] + tmp(1,1)*v[NS*ix+1];
	v.push_back( tmp2(1,0) ); rows.push_back( NS*ix+1 ); cols.push_back( NS*ix );
	v.push_back( tmp2(1,1) ); rows.push_back( NS*ix+1 ); cols.push_back( NS*ix+1 );
      }
    }
  }



template <typename T>
void matmul( T* res, T* v,
	     const std::vector<Complex>& val,
	     const std::vector<int>& cols,
	     const std::vector<int>& rows ) {
  const int N = rows.size();

  for(int i=0; i<N; i++){
    res[i] = 0.0;
    const int row_start = rows[i];
    const int row_end = rows[i+1];
    for(int jj=row_start; jj<row_end; jj++){
      res[i] += val[jj] * v[cols[jj]];
    }
  }
}


template <typename T>
void matmulcoo( T* res, T* v,
		const std::vector<Complex>& val,
		const std::vector<int>& cols,
		const std::vector<int>& rows,
		const int N) {
  assert( val.size()==cols.size() );
  assert( val.size()==rows.size() );

  for(int i=0; i<N; i++) res[i] = 0.0;
  for(int i=0; i<val.size(); i++) res[rows[i]] += val[i] * v[cols[i]];
}

template <typename T=Complex>
void matmuladjointcoo( T* res, T* v,
		       const std::vector<Complex>& val,
		       const std::vector<int>& cols,
		       const std::vector<int>& rows,
		       const int N) {
  assert( val.size()==cols.size() );
  assert( val.size()==rows.size() );

  for(int i=0; i<N; i++) res[i] = 0.0;
  for(int i=0; i<val.size(); i++) res[cols[i]] += std::conj(val[i]) * v[rows[i]];
}


template <typename T>
void matmulgam5( T* res, T* v, const int Nx) {
  for(int ix=0; ix<Nx; ix++){
    res[2*ix] = v[2*ix];
    res[2*ix+1] = -v[2*ix+1];
  }
}


  {
    U1onS2 ds = SW.d(U);

    const double eps = 1.0e-2;

    for(int a=0; a<U.field.size(); a++){
      U1onS2 UP(U);
      U1onS2 UM(U);

      UP[a] += eps;
      UM[a] -= eps;

      double check = ( SW(UP) - SW(UM) ) / (2.0*eps);
      std::cout << "ds = " << ds[a] << std::endl
                << "ch = " << check << std::endl
		<< "ratio = " << ds[a]/check << std::endl
		<< std::endl;
    }
  }




  // {
  //   auto mat = D.matrix_form( U );
  //   // auto gam5_D = matmultgam5( mat );
  //   // std::cout << gam5_D.adjoint() - gam5_D << std::endl;
  //   // std::cout << "det D = "  << mat.determinant() << std::endl;
  //   // std::cout << "det HW = " << gam5_D.determinant() << std::endl;

  //   // std::cout << "det = " << mat.determinant() << std::endl;

  //   // // ----------------

  //   const CGCUDA cg( D );

  //   Complex v[cg.sparse.N];
  //   for(int i=0; i<cg.sparse.N; i++) v[i] = rng.gaussian();

  //   Complex res[cg.sparse.N];
  //   cg( res, v, U );

  //   const VC r = Eigen::Map<VC>( v, cg.sparse.N );
  //   auto tmp1 = (mat.adjoint()*mat).inverse() * r;

  //   double norm = 0.0;
  //   for(int i=0; i<cg.sparse.N; i++) {
  //     // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
  //     norm += std::abs(tmp1[i] - res[i]);
  //   }
  //   std::cout << "norm = " << norm << std::endl;
  // }


  // {
  //   const CGCUDA cg( D );

  //   Complex v[cg.sparse.N];
  //   std::vector<Complex> xi(cg.sparse.N, 0.0);
  //   std::vector<Complex> phi(cg.sparse.N, 0.0);
  //   for(int i=0; i<cg.sparse.N; i++) xi[i] = rng.gaussian();

  //   // Complex res[cg.sparse.N];
  //   Complex D_coo[cg.sparse.len], D_csrH[cg.sparse.len];
  //   D.coo_format(D_coo, U);
  //   cg.sparse.coo2csrH( D_csrH, D_coo );
  //   cg.sparse.multT<Complex>( phi.data(), xi.data(), D_csrH );

  //   const VC r = Eigen::Map<VC>( xi.data(), cg.sparse.N );
  //   auto mat = D.matrix_form( U );
  //   auto tmp1 = mat.adjoint() * r;
  //   for(int i=0; i<phi.size(); i++) std::cout << "diff = " << tmp1(i) - phi[i] << std::endl;
  // }


  // for(int ix=0; ix<D.lattice.n_sites; ix++){
  //   std::cout << ix << " " << rng.gaussian_site(ix) << std::endl;
  // }


  {
    PseudoFermion phi( D );
    phi.gen( U, rng );

    Eigen::MatrixXcd tmp = Eigen::MatrixXcd::Zero( phi.phi.size(), phi.phi.size() );
    int kkmax = 1000000;

    for(int kk=0; kk<kkmax; kk++){
      phi.gen( U, rng );
      for(int i=0; i<phi.phi.size(); i++) for(int j=0; j<phi.phi.size(); j++) {
	  // tmp(i,j) += phi[i] * std::conj(phi[j]);
	  // tmp(i,j) += std::conj(phi[i]) * phi[j];
	  tmp(i,j) += phi[i] * std::conj(phi[j]);
	  // std::cout << *iter << " ";
	}
    }
    tmp /= kkmax;

    // for(int i=0; i<phi.phi.size(); i++) {
    //   for(int j=0; j<phi.phi.size(); j++) {
    // 	std::cout << i << " " << j << " " << tmp(i,j) << " ";
    //   }
    //   std::cout << std::endl;
    // }
    // std::cout << std::endl;

    auto mat = D.matrix_form( U );
    auto mat2= mat.adjoint() * mat;
    std::cout << "diff (Ddag D) = " << std::endl;
    // auto mat2= mat * mat.adjoint();
    // std::cout << "D Ddag = " << std::endl;
    for(int i=0; i<phi.phi.size(); i++) {
      for(int j=0; j<phi.phi.size(); j++) {
	std::cout << i << " " << j << " " << mat2(i,j) << " " << tmp(i,j) << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;


    // auto mat = D.matrix_form( U );
    // auto mat2= mat.adjoint() * mat;
    // std::cout << "diff (Ddag D) = " << std::endl;
    // // auto mat2= mat * mat.adjoint();
    // // std::cout << "D Ddag = " << std::endl;
    // for(int i=0; i<phi.phi.size(); i++) {
    //   for(int j=0; j<phi.phi.size(); j++) {
    // 	std::cout << i << " " << j << " " << mat2(i,j) - tmp(j,i) << std::endl;
    //   }
    //   std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout << mat2 << std::endl;
    // std::cout << mat2.inverse() << std::endl;
  }






  PseudoFermion pf( D );
  pf.gen( U, rng );

  {
    using Link = std::array<int,2>; // <int,int>;

    std::vector<Complex> dD;
    std::vector<int> is;
    std::vector<int> js;

    const int ix=0, jj=0;
    const int iy=lattice.nns[ix][jj];

    D.d_coo_format( dD, is, js, U, Link{ix,iy} );

    // std::cout << "vD = " << std::endl;
    // for(int i=0; i<vD.size(); i++) {
    //   std::cout << is[i] << " " << js[i] << " " << vD[i] << std::endl;
    // }
    // std::cout << std::endl;
    // std::cout << "is = " << std::endl;
    // for(int i=0; i<is.size(); i++) std::cout << ;
    // std::cout << std::endl;
    // std::cout << "js = " << std::endl;
    // for(int i=0; i<js.size(); i++) std::cout << js[i] << " ";
    // std::cout << std::endl;

    std::vector<Complex> v(lattice.n_sites*2, 0.0);
    constexpr Complex I = Complex(0.0, 1.0);
    for(int ix=0; ix<lattice.n_sites; ix++) for(int a=0; a<2; a++) v[2*ix+a] = ( rng.gaussian_site(ix) + I*rng.gaussian_site(ix) ) / std::sqrt(2.0);

    std::vector<Complex> dDv;
    pf.cg.sparse.multcoo( dDv, v, dD, is, js );

    std::cout << "dDv = " << std::endl;
    for(int i=0; i<dDv.size(); i++) std::cout << dDv[i] << " ";
    std::cout << std::endl;

    const double eps = 1.0e-5;

    // for(int a=0; a<U.field.size(); a++){
    const int ell=lattice.map2il.at(Link{ix,iy});
    Gauge UP(U);
    Gauge UM(U);

    UP[ell] += eps;
    UM[ell] -= eps;

    auto DP = D.matrix_form( UP );
    auto DM = D.matrix_form( UM );

    const VC r = Eigen::Map<VC>( v.data(), pf.cg.sparse.N );
    // auto numeric = ( DP - DM ) / (2.0*eps);
    // for(int i=0; i<numeric.cols(); i++) {
    //   for(int j=0; j<numeric.rows(); j++) {
    // 	if( std::abs(numeric(i,j))>1.0e-15 ) std::cout << i << " " << j << " " << numeric(i,j) << std::endl;
    //   }
    //   // std::cout << std::endl;
    // }
    // std::cout << std::endl;

    auto numeric = ( DP*r - DM*r ) / (2.0*eps);
    std::cout << "numeric = " << std::endl;
    for(int i=0; i<numeric.size(); i++) std::cout << numeric[i] << " ";
    std::cout << std::endl;
  }





  {
    PseudoFermion phi( D );

    Eigen::MatrixXcd tmp = Eigen::MatrixXcd::Zero( phi.phi.size(), phi.phi.size() );
    int kkmax = 1000;

    for(int kk=0; kk<kkmax; kk++){
      phi.gen( U, rng );
      std::vector<Complex> eta = phi.get_eta( U );
      for(int i=0; i<phi.phi.size(); i++) for(int j=0; j<phi.phi.size(); j++) {
	  tmp(i,j) += eta[i] * std::conj(eta[j]);
	}
    }
    tmp /= kkmax;

    auto mat = D.matrix_form( U );
    auto mat2= (mat.adjoint() * mat).inverse();
    std::cout << "diff (Ddag D)^{-1} = " << std::endl;
    for(int i=0; i<phi.phi.size(); i++) {
      for(int j=0; j<phi.phi.size(); j++) {
	std::cout << i << " " << j << " " << mat2(i,j) << " " << tmp(i,j) << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }



  {
    PseudoFermion phi( D );

    Eigen::MatrixXcd tmp = Eigen::MatrixXcd::Zero( phi.phi.size(), phi.phi.size() );
    int kkmax = 1000;

    for(int kk=0; kk<kkmax; kk++){
      phi.gen( U, rng );
      std::vector<Complex> eta = phi.get_eta( U );
      for(int i=0; i<phi.phi.size(); i++) for(int j=0; j<phi.phi.size(); j++) {
	  tmp(i,j) += eta[i] * std::conj(eta[j]);
	}
    }
    tmp /= kkmax;

    auto mat = D.matrix_form( U );
    auto mat2= (mat.adjoint() * mat).inverse();
    std::cout << "diff (Ddag D)^{-1} = " << std::endl;
    for(int i=0; i<phi.phi.size(); i++) {
      for(int j=0; j<phi.phi.size(); j++) {
	std::cout << i << " " << j << " " << mat2(i,j) << " " << tmp(i,j) << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }


// VE circumcenter(const VE& r0, const VE& r1, const VE& r2){
//   const VE r10 = r1 - r0;
//   const VE r20 = r2 - r0;

//   const VE tmp1 = r10.squaredNorm() * r20 - r20.squaredNorm() * r10;
//   const VE cross = r10.cross(r20);
//   const VE numer = tmp1.cross(cross);
//   const double denom = 2.0*cross.squaredNorm();

//   return numer/denom + r0;
// }


  // MS gamma(const int ix, const int iy, const double shift=0.0) const { // located at x
  //   const double al = alpha.at(Link{ix,iy}) + shift;
  //   // return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  //   return std::cos(al)*sigma[1] + std::sin(al)*sigma[2];
  // }





  // VC operator()( const VC& v ) const {
  //   VC res = VC::Zero(v.size());

  //   for(int ix=0; ix<lattice.n_sites; ix++){
  //     for(int jj=0; jj<3; jj++){
  // 	int iy = lattice.nns[ix][jj];
	
  // 	{
  // 	  // res.block<NS,NS>(NS*ix,NS*iy) -= lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * Omega(ix, iy);
  // 	  const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * Omega(ix, iy);
  // 	  res.segment<NS>(NS*ix) -= tmp * v.segment<NS>(NS*iy);
  // 	}

  // 	{
  // 	  // res.block<NS,NS>(NS*ix,NS*ix) += lattice.vol[ix]/lattice.mean_vol * r*lattice.u[ix][jj]*sigma[0];
  // 	  const MS tmp = lattice.vol[ix]/lattice.mean_vol * r*lattice.u[ix][jj]*sigma[0];
  // 	  res.segment<NS>(NS*ix) += tmp * v.segment<NS>(NS*ix);
  // 	}
  //     }
  //   }

  //   return res;
  // } // end matrix_form


  // VC operator()( const U1onS2& U, const VC& v ) const {
  //   VC res = VC::Zero(v.size());

  //   for(int ix=0; ix<lattice.n_sites; ix++){
  //     for(int jj=0; jj<3; jj++){
  // 	int iy = lattice.nns[ix][jj];
	
  // 	{
  // 	  const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);
  // 	  res.segment<NS>(NS*ix) -= tmp * v.segment<NS>(NS*iy);
  // 	}

  // 	{
  // 	  const MS tmp = lattice.vol[ix]/lattice.mean_vol * r*lattice.u[ix][jj]*sigma[0];
  // 	  res.segment<NS>(NS*ix) += tmp * v.segment<NS>(NS*ix);
  // 	}
  //     }
  //   }

  //   return res;
  // } // end matrix_form





  {
    using Link = std::array<int,2>; // <int,int>;
    PseudoFermion phi( D, U, rng );
    Force exact = phi.get_force( U );

    for(int il=0; il<U.lattice.n_links; il++){
      const double eps = 1.0e-5;
      Gauge UP(U);
      Gauge UM(U);

      UP[il] += eps;
      UM[il] -= eps;

      std::vector<Complex> etaP = phi.get_eta( UP );
      std::vector<Complex> etaM = phi.get_eta( UM );

      Complex SP = phi.dot( phi.phi, etaP );
      Complex SM = phi.dot( phi.phi, etaM );
      Complex numeric = ( SP - SM ) / (2.0*eps);
      std::cout << exact[il] << " " << numeric << std::endl;
    }
  }



  {
    using Link = std::array<int,2>; // <int,int>;
    PseudoFermion phi( D, U, rng );
    Force exact = phi.get_force( U );

    for(int il=0; il<U.lattice.n_links; il++){
      const double eps = 1.0e-5;
      Gauge UP(U);
      Gauge UM(U);

      UP[il] += eps;
      UM[il] -= eps;

      double numeric = ( phi.S(UP) - phi.S(UM) ) / (2.0*eps);
      std::cout << exact[il] << " " << numeric << std::endl;
    }
  }



  {
    using Link = std::array<int,2>; // <int,int>;
    PseudoFermion phi( D, U, rng );
    Force exact = phi.dS( U );

    for(int il=0; il<U.lattice.n_links; il++){
      const double eps = 1.0e-5;
      Gauge UP(U);
      Gauge UM(U);

      UP[il] += eps;
      UM[il] -= eps;

      double numeric = ( phi.S(UP) - phi.S(UM) ) / (2.0*eps);
      std::cout << exact[il] << " " << numeric << std::endl;
    }
  }



  double stot = 1.0;
  // int nsteps = 10;

  // ---------------------------------------

  Force pi( lattice );
  pi.gaussian( rng );

  for(int nsteps=10; nsteps<100; nsteps+=10){
    HMC hmc(rng, SW, D, stot, nsteps);
    // HMC<GaugeForce,GaugeField,GaugeAction> hmc(rng, SW, stot, nsteps);

    Gauge U1( U );
    rng.reseed( 1 );
    PseudoFermion phi( D, U, rng );
    Force pi1(pi);

    // std::cout << "pi1 = " << std::endl;
    // for(auto elem : pi1 ) std::cout << elem << " ";
    // std::cout << std::endl;
    // std::cout << "U1 = " << std::endl;
    // for(auto elem : U1 ) std::cout << elem << " ";
    // std::cout << std::endl;

    const double h0 = hmc.H(pi1, U1, phi);
    hmc.leapfrog_explicit( pi1, U1, phi );
    const double h1 = hmc.H(pi1, U1, phi);

    std::cout << nsteps << " " << h1-h0 << std::endl;
  }  

