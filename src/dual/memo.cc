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




  // double stot = 1.0;

  // Force pi( lattice );
  // pi.gaussian( rng );

  // for(int nsteps=10; nsteps<100; nsteps+=10){
  //   rng.reseed( 1 );
  //   HMC hmc(rng, SW, D, stot, nsteps);

  //   Force pi1(pi);
  //   Gauge U1( U );
  //   hmc.phi.gen( U, rng );

  //   const double h0 = hmc.H(pi1, U1);
  //   hmc.leapfrog_explicit( pi1, U1 );
  //   const double h1 = hmc.H(pi1, U1);

  //   std::cout << nsteps << " " << h1-h0 << std::endl;
  // }

  // ---------------------------------------
<<<<<<< HEAD





  // {
  //   double z = 0.1;

  //   double k, sn, cn, K;
  //   z = 0.1;
  //   k = 0.2;
  //   EllipticUtils::sncnK( z, k, sn, cn, K);
  //   std::cout << "z  = " << z  << std::endl
  // 	      << "k  = " << k  << std::endl
  // 	      << "sn = " << sn << std::endl
  // 	      << "cn = " << cn << std::endl
  // 	      << "K  = " << K  << std::endl;

  //   double kp, snk, cnk, Kp;
  //   kp = std::sqrt(1.0-k*k);
  //   EllipticUtils::sncnK( z, kp, snk, cnk, Kp, agm );
  //   std::cout << "Kp = " << Kp << std::endl;
  // }


  {
    const int n=11;
    const double k = 0.1;

    const double kp = std::sqrt(1.0-k*k);

    std::vector<double> c( int(n/2)+1, 0.0 );
    std::vector<double> cp( int(n/2)+1, 0.0 );

    double Kp = 1.0, xibar;
    for(int m=0; m<=int(n/2); m++){
      double sn, cn, dn;
      double z = 2.0 * Kp * m / n;
      EllipticUtils::sncndnK( z, kp, sn, cn, dn, Kp );
      if(m==0) continue;
      c[m] = - std::pow( cn / sn, 2 );

      z = 2.0 * Kp * (m-0.5) / n;
      EllipticUtils::sncndnK( z, kp, sn, cn, dn, Kp );
      cp[m] = - std::pow( cn / sn, 2 );
      if(m==1) xibar = 1.0/dn;

      std::cout << "c[" << m << "] = " << c[m] << std::endl
		<< "cp[" << m << "] = " << cp[m] << std::endl;
    }

    double M = 1.0;
    for(int m=1; m<=int(n/2); m++) M *= (1.0-c[m]) / (1.0-cp[m]);

    double lambda_inv = xibar / M;
    for(int m=1; m<=int(n/2); m++) lambda_inv *= (1.0-c[m]*xibar*xibar) / (1.0-cp[m]*xibar*xibar);

    // double A = 2.0 / (1.0+1.0/lambda) / k;
    // for(int m=1; m<=int(n/2); m++) A *= ap[m] / a[m] * (1.0-k*k/ap[m]) / (1.0-k*k/a[m]);
    std::cout << "M = " << M << std::endl
	      << "lambda_inv = " << lambda_inv << std::endl;

    {
      for(double x = 0.01; x<=1.0; x+=0.01){
	double res = 2.0 / (1.0+lambda_inv) * x / (k*M);
	for(int m=1; m<=int(n/2); m++) res *= (k*k - c[m]*x*x) / (k*k - cp[m]*x*x);
	std::cout << x << " " << 1.0-res << std::endl;
      }
    }

    const double Delta = (lambda_inv - 1.0) / (lambda_inv + 1.0);
    std::cout << "Delta = " << Delta << std::endl;
  }



  Zolotarev f;
  std::cout << "Delta = " << f.Delta() << std::endl;

  for(double x = -1.0; x<=1.0; x+=0.01){
    std::cout << x << " " << 1.0-std::abs(f(x)) << std::endl;
  }

// struct CGCUDA{
//   const int N;
//   const int len;

//   CuC *d_val, *d_valH;
//   int *d_cols, *d_rows, *d_colsT, *d_rowsT;
//   CuC *d_x, *d_r, *d_p, *d_q, *d_tmp;
//   CuC *d_scalar;


//   CGCUDA( const Sparse& sparse )
//     : N(sparse.N)
//     , len(sparse.len)
//   {
//     int device;
//     cudacheck(cudaGetDeviceCount(&device));
//     cudaDeviceProp device_prop[device];
//     cudaGetDeviceProperties(&device_prop[0], 0);
//     cudacheck(cudaSetDevice(0));// "TITAN V"

//     // sparse matrix
//     cudacheck(cudaMalloc(&d_val, len*CD));
//     cudacheck(cudaMalloc(&d_valH, len*CD));
//     //
//     cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
//     cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
//     cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
//     cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));

//     cudacheck(cudaMemcpy(d_cols, sparse.cols_csr.data(), len*sizeof(int), H2D));
//     cudacheck(cudaMemcpy(d_rows, sparse.rows_csr.data(), (N+1)*sizeof(int), H2D));
//     cudacheck(cudaMemcpy(d_colsT, sparse.cols_csrT.data(), len*sizeof(int), H2D));
//     cudacheck(cudaMemcpy(d_rowsT, sparse.rows_csrT.data(), (N+1)*sizeof(int), H2D));

//     // ------------------

//     // CG
//     cudacheck(cudaMalloc(&d_x, N*CD));
//     cudacheck(cudaMalloc(&d_r, N*CD));
//     cudacheck(cudaMalloc(&d_p, N*CD));
//     cudacheck(cudaMalloc(&d_q, N*CD));
//     cudacheck(cudaMalloc(&d_tmp, N*CD));

//     cudacheck(cudaMalloc(&d_scalar, CD));

//   }

//   ~CGCUDA(){
//     cudacheck(cudaFree(d_x));
//     cudacheck(cudaFree(d_r));
//     cudacheck(cudaFree(d_p));
//     cudacheck(cudaFree(d_q));
//     cudacheck(cudaFree(d_tmp));
//     cudacheck(cudaFree(d_scalar));

//     cudacheck(cudaFree(d_val));
//     cudacheck(cudaFree(d_valH));
//     cudacheck(cudaFree(d_cols));
//     cudacheck(cudaFree(d_rows));
//     cudacheck(cudaFree(d_colsT));
//     cudacheck(cudaFree(d_rowsT));

//     cudacheck(cudaDeviceReset());
//   }
  


//   __host__
//   void solve(CuC* x, CuC* b,
// 	     CuC* val, CuC* valH,
// 	     const double tol=1.0e-13, const int maxiter=1e8){
//     // CG
//     cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_r, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_p, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_q, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_tmp, 0, N*CD)); // added @@
//     cudacheck(cudaMemset(d_scalar, 0, CD)); // added @@

//     cudacheck(cudaMemcpy(d_r, b, N*CD, H2D));
//     cudacheck(cudaMemcpy(d_p, d_r, N*CD, D2D));

//     double mu; dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
//     assert(mu>=0.0);
//     double mu_old = mu;

//     double b_norm_sq; dot2self_normalized_wrapper(b_norm_sq, d_scalar, d_r, N);
//     assert(b_norm_sq>=0.0);
//     double mu_crit = tol*tol*b_norm_sq;

//     if(mu<mu_crit) std::clog << "NO SOLVE" << std::endl;
//     else{
//       int k=0;
//       CuC gam;

//       for(; k<maxiter; ++k){
// 	multA(d_q, d_tmp, d_p,
// 	      d_val, d_cols, d_rows,
// 	      d_valH, d_colsT, d_rowsT,
// 	      N
// 	      );

// 	dot_normalized_wrapper(gam, d_scalar, d_p, d_q, N);

// 	CuC al = mu/gam;
// 	cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
// 	daxpy<<<NBlocks, NThreadsPerBlock>>>(d_x, d_scalar, d_p, d_x, N);

// 	al = -al;
// 	cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
// 	daxpy<<<NBlocks, NThreadsPerBlock>>>(d_r, d_scalar, d_q, d_r, N);

// 	dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
// 	assert(mu>=0.0);

// 	if(mu<mu_crit || std::isnan(mu)) break;
// 	CuC bet = cplx(mu/mu_old);
// 	mu_old = mu;

// 	cudacheck(cudaMemcpy(d_scalar, &bet, CD, H2D));
// 	daxpy<<<NBlocks, NThreadsPerBlock>>>(d_p, d_scalar, d_p, d_r, N);

// 	if(k%100==0) {
// 	  std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
// 	}
//       }
//       std::clog << "SOLVER:       #iterations: " << k << std::endl;
//       std::clog << "SOLVER:       mu =         " << mu << std::endl;
//     }

//     cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));
//   }

// };


// __global__
// void mult_gam5_M( CuC* res,
// 		  const CuC* v,
// 		  const CuC* v_csr,
// 		  const int* cols,
// 		  const int* rows,
// 		  const int N
// 		  ){
//   Idx i = blockIdx.x*blockDim.x + threadIdx.x;

//   if(i<N) {
//     res[i] = cplx(0.0);
//     const int row_start = rows[i];
//     const int row_end = rows[i+1];
//     for(int jj=row_start; jj<row_end; jj++) res[i] = res[i] + v_csr[jj] * v[ cols[jj] ];
//     res[i] = ( -2*(i%2) + 1 ) * res[i];
//   }
// }


// __global__
// void multgam5_self ( CuC* v,
// 		     const int N
// 		){
//   Idx i = blockIdx.x*blockDim.x + threadIdx.x;
//   if(i<N) v[i] = ( -2*(i%2) + 1 ) * v[i];
// }


// __host__
// void multA(CuC* d_v, CuC* d_tmp, CuC* d_v0,
// 	   CuC* d_val, int* d_cols, int* d_rows,
// 	   // CuC* d_valH, int* d_colsT, int* d_rowsT,
// 	   const int N
// 	   ){
//   cudacheck(cudaMemset(d_tmp, 0, N*CD));
//   mult<<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v0, d_val, d_cols, d_rows, N);

//   cudacheck(cudaMemset(d_v, 0, N*CD));
//   // mult_gam5_M<<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v0, d_val, d_cols, d_rows, N);
//   mult<<<NBlocks, NThreadsPerBlock>>>(d_v, d_tmp, d_valH, d_colsT, d_rowsT, N);
// }



// __host__
// void solve(CuC* x, const CuC* b,
// 	   const CuC* val, const std::vector<int>& cols, const std::vector<int>& rows,
// 	   // const CuC* valH, const std::vector<int>& colsT, const std::vector<int>& rowsT,
// 	   const int N, const int len,
// 	   const double tol=1.0e-13, const int maxiter=1e8){

//   // sparse matrix
//   CuC *d_val, *d_valH;
//   cudacheck(cudaMalloc(&d_val, len*CD));
//   cudacheck(cudaMemcpy(d_val, val, len*CD, H2D));
//   //
//   // cudacheck(cudaMalloc(&d_valH, len*CD));
//   // cudacheck(cudaMemcpy(d_valH, valH, len*CD, H2D));
//   //
//   int *d_cols, *d_rows; // , *d_colsT, *d_rowsT;
//   cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
//   cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
//   cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(int), H2D));
//   cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(int), H2D));
//   //
//   // cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
//   // cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));
//   // cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(int), H2D));
//   // cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(int), H2D));

//   // CG
//   CuC *d_x, *d_r, *d_p, *d_q, *d_tmp;
//   cudacheck(cudaMalloc(&d_x, N*CD));
//   cudacheck(cudaMalloc(&d_r, N*CD));
//   cudacheck(cudaMalloc(&d_p, N*CD));
//   cudacheck(cudaMalloc(&d_q, N*CD));
//   cudacheck(cudaMalloc(&d_tmp, N*CD));
//   cudacheck(cudaMemset(d_x, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_r, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_p, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_q, 0, N*CD)); // added @@
//   cudacheck(cudaMemset(d_tmp, 0, N*CD)); // added @@

//   CuC *d_scalar;
//   cudacheck(cudaMalloc(&d_scalar, CD));
//   cudacheck(cudaMemset(d_scalar, 0, CD)); // added @@

//   cudacheck(cudaMemcpy(d_r, b, N*CD, H2D));
//   cudacheck(cudaMemcpy(d_p, d_r, N*CD, D2D));

//   double mu; dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
//   assert(mu>=0.0);
//   double mu_old = mu;

//   double b_norm_sq; dot2self_normalized_wrapper(b_norm_sq, d_scalar, d_r, N);
//   assert(b_norm_sq>=0.0);
//   double mu_crit = tol*tol*b_norm_sq;

//   if(mu<mu_crit) std::clog << "NO SOLVE" << std::endl;
//   else{
//     int k=0;
//     CuC gam;

//     for(; k<maxiter; ++k){
//       // multA(d_q, d_tmp, d_p, nu);
//       multA(d_q, d_tmp, d_p,
// 	    d_val, d_cols, d_rows,
// 	    // d_valH, d_colsT, d_rowsT,
// 	    N
// 	    );

//       dot_normalized_wrapper(gam, d_scalar, d_p, d_q, N);

//       CuC al = mu/gam;
//       cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
//       daxpy<<<NBlocks, NThreadsPerBlock>>>(d_x, d_scalar, d_p, d_x, N);

//       al = -al;
//       cudacheck(cudaMemcpy(d_scalar, &al, CD, H2D));
//       daxpy<<<NBlocks, NThreadsPerBlock>>>(d_r, d_scalar, d_q, d_r, N);

//       dot2self_normalized_wrapper(mu, d_scalar, d_r, N);
//       assert(mu>=0.0);

//       if(mu<mu_crit || std::isnan(mu)) break;
//       CuC bet = cplx(mu/mu_old);
//       mu_old = mu;

//       cudacheck(cudaMemcpy(d_scalar, &bet, CD, H2D));
//       daxpy<<<NBlocks, NThreadsPerBlock>>>(d_p, d_scalar, d_p, d_r, N);

//       if(k%100==0) {
// 	std::clog << "SOLVER:       #iterations: " << k << ", mu =         " << mu << std::endl;
//       }
//     }
//     std::clog << "SOLVER:       #iterations: " << k << std::endl;
//     std::clog << "SOLVER:       mu =         " << mu << std::endl;
//   }

//   cudacheck(cudaMemcpy(x, d_x, N*CD, D2H));

//   cudacheck(cudaFree(d_x));
//   cudacheck(cudaFree(d_r));
//   cudacheck(cudaFree(d_p));
//   cudacheck(cudaFree(d_q));
//   cudacheck(cudaFree(d_tmp));
//   cudacheck(cudaFree(d_scalar));

//   cudacheck(cudaFree(d_val));
//   cudacheck(cudaFree(d_cols));
//   cudacheck(cudaFree(d_rows));
//   // cudacheck(cudaFree(d_valH));
//   // cudacheck(cudaFree(d_colsT));
//   // cudacheck(cudaFree(d_rowsT));
// }



    


  // gpu
  SparseHelper H_DW(lattice, true);
  H_DW.set_dirac();
  SparseMatrix<CuC> M_DW, M_DWH;
  H_DW.reset_DU( DW, U );
  M_DW.assign_from_helper_DW_gpu(H_DW);
  M_DWH.assign_from_helper_DWH_gpu(H_DW);

  // LinOp<CuC> Op_DW;
  // Op_DW.set_coeff ( 0,    cplx(1.0) );
  // Op_DW.set_matrix( 0, 0, &M_DW );

  // LinOp<CuC> Op_DWHDW({2});
  // Op_DWHDW.set_coeff ( 0,    cplx(1.0) );
  // Op_DWHDW.set_matrix( 0, 0, &M_DW );
  // Op_DWHDW.set_matrix( 0, 1, &M_DWH );

  LinOp<CuC> Op_DWHDW_const({2, 0});
  const CuC a = cplx(-2.0);
  Op_DWHDW_const.set_coeff ( 0,    cplx(1.0) );
  Op_DWHDW_const.set_matrix( 0, 0, &M_DW );
  Op_DWHDW_const.set_matrix( 0, 1, &M_DWH );
  Op_DWHDW_const.set_coeff ( 1,    a );

  constexpr Idx N = CompilationConst::N;
  Eigen::MatrixXcd mat(N, N);
  {
    // auto mat = Dov.D.matrix_form();
    // auto mat = DW.matrix_form();

    double lambda_max = 6.0;

    for(Idx i=0; i<N; i++){
      Eigen::VectorXcd e = Eigen::VectorXcd::Zero(N);
      e(i) = 1.0;
      std::vector<Complex> xi(e.data(), e.data()+N);
      std::vector<Complex> Dxi(N);

      // H_DW.reset_DU( DW, U );
      // Op_DW.from_cpu<N>( Dxi, xi );
      // Op_DWHDW.from_cpu<N>( Dxi, xi );
      Op_DWHDW_const.from_cpu<N>( Dxi, xi );
      // Op_DW.on_cpu<N>( Dxi, xi );
      
      // Dov.multHW( Dxi, xi, U, lambda_max );
      // matmulgam5<Complex>( Dxi.data(), Dxi.data(), int(N/2) );
      // mult_a<Complex>( Dxi.data(), 6.0, N );

      // Dov.H( Dxi.data(), xi.data(), U, 6.0 );
      // mult_a<Complex>( Dxi.data(), 6.0, N );

      mat.block(0,i,N,1) = Eigen::Map<Eigen::MatrixXcd>(Dxi.data(), N, 1);
    }
  }


  // void set_diag(){
  //   len = N;

  //   // ========= CSR ========= //

  //   cols_csr.resize(N);
  //   for(Idx i=0; i<N; i++){
  //     cols_csr[i] = i;
  //     rows_csr.push_back( i );
  //   }
  //   rows_csr.push_back( N );

  //   cols_csrT = cols_csr;
  //   rows_csrT = rows_csr;

  //   assert( rows_csr.size()==N+1 );
  //   assert( rows_csrT.size()==N+1 );

  //   if(locate_on_gpu){
  //     const std::vector<Idx>& cols = cols_csr;
  //     const std::vector<Idx>& rows = rows_csr;
  //     const std::vector<Idx>& colsT = cols_csrT;
  //     const std::vector<Idx>& rowsT = rows_csrT;

  //     cudacheck(cudaMalloc(&d_cols, len*sizeof(Idx)));
  //     cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
  //     cudacheck(cudaMemcpy(d_cols, cols.data(), len*sizeof(Idx), H2D));
  //     cudacheck(cudaMemcpy(d_rows, rows.data(), (N+1)*sizeof(Idx), H2D));
  //     //
  //     cudacheck(cudaMalloc(&d_colsT, len*sizeof(Idx)));
  //     cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(Idx)));
  //     cudacheck(cudaMemcpy(d_colsT, colsT.data(), len*sizeof(Idx), H2D));
  //     cudacheck(cudaMemcpy(d_rowsT, rowsT.data(), (N+1)*sizeof(Idx), H2D));
  //   }

  //   v_coo.resize(len);
  //   v_csr.resize(len);
  //   v_csrH.resize(len);

  //   is_set = true;
  // }
  const double M5 = -2.0;
  WilsonDirac DW(lattice, M5);

  SparseMatrix<CuC> M_DW, M_DWH;

  // double lambda_max = 10.0;
  // MatPoly<CuC> Op;
  // // Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  // const CuC a = cplx(-0.5);
  // Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  // Op.push_back ( a/(lambda_max*lambda_max), {} );


  Overlap Dov;


  constexpr Idx N = CompilationConst::N;
  Eigen::MatrixXcd mat(N, N);
  {
    for(Idx i=0; i<N; i++){
      Eigen::VectorXcd e = Eigen::VectorXcd::Zero(N);
      e(i) = 1.0;
      std::vector<Complex> xi(e.data(), e.data()+N);
      std::vector<Complex> Dxi(N);

      // Op.solve<N>( Dxi, xi );
      // xi = Dxi;
      // Op.from_cpu<N>( Dxi, xi );
      Dov( Dxi, xi, &M_DW, &M_DWH );

      mat.block(0,i,N,1) = Eigen::Map<Eigen::MatrixXcd>(Dxi.data(), N, 1);
    }
  }




  std::vector<Complex> q(N, 0.0);
  std::vector<Complex> x(N, 0.0);
  // std::vector<Complex> eta(N, 0.0);
  // xi[0] = 0.4;
  // xi[1] = -0.2;
  // #ifdef _OPENMP
  // #pragma omp parallel for num_threads(CompilationConst::NPARALLEL)
  // #endif
  for(int i=0; i<N; i++) {
    std::cout << "i = " << i << std::endl;
    // xi[i] = rng.gaussian_site(i);
    q[i] = rng.gaussian();
    // x[i] = rng.gaussian();
  }

  std::cout << "init set." << std::endl;

  {
    MatPoly<CuC> Op;
    Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );

    CuC *d_x; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_x, N*CD));
    CuC *d_q; // , *d_tmp, *d_tmp2;
    CUDA_CHECK(cudaMalloc(&d_q, N*CD));
    CuC *d_scalar;
    CUDA_CHECK(cudaMalloc(&d_scalar, CD));

    Complex dot;
    double norm=1.0, mu_0=1.0, mu_m1=1.0, mu_m2=1.0;

    CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));
    Op.dot2self<N>(norm, d_scalar, d_q);
    for(int i=0; i<N; i++) q[i] = q[i]/std::sqrt(norm);

    double TOL=1.0e-4;

    double lambda=100.0, lambda_old=1000.0;

    for(int i=0; i<200; i++){
      std::vector<Complex> DHDxi(N);

      Op.from_cpu<N>( x, q );
      CUDA_CHECK(cudaMemcpy(d_x, reinterpret_cast<const CuC*>(x.data()), N*CD, H2D));

      //

      Op.dot2self<N>(norm, d_scalar, d_x);
      for(int i=0; i<N; i++) q[i] = x[i]/std::sqrt(norm);
      CUDA_CHECK(cudaMemcpy(d_q, reinterpret_cast<const CuC*>(q.data()), N*CD, H2D));
      
      Op.dot<N>(reinterpret_cast<CuC&>(dot), d_scalar, d_x, d_q);
      // CUDA_CHECK(cudaMemcpy(d_Dxi, reinterpret_cast<const CuC*>(DHDxi.data()), N*CD, H2D));

      // std::cout << std::sqrt(norm) << ", " << DHDxi[0] << std::endl;
      mu_m2=mu_m1;
      mu_m1=mu_0;
      mu_0=dot.real();
      // std::cout << dot.real() << std::endl;
      double r = (mu_0-mu_m1)/(mu_m1-mu_m2);
      double a = (mu_0-mu_m1)/std::pow(r,i-1)/(r-1);
      lambda_old = lambda;
      lambda = mu_0 - a*std::pow(r,i);
      std::cout << i << " " << dot.real() << " " << lambda << std::endl;

      if(std::abs(lambda_old-lambda)/lambda<TOL) break;
      // std::cout << dot.real() << std::endl;
      // dot_old = dot;

      // for(int i=0; i<N; i++) xi[i] = DHDxi[i];
      // for(int i=0; i<N; i++) xi[i] = DHDxi[i]/dot;
      // double norm;
      // CUDA_CHECK(cudaMemcpy(d_Dxi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));
      // Op.dot2self<N>(norm, d_scalar, d_Dxi);
      //std::cout << norm << ", " << DHDxi[0] << std::endl;
    }

    CUDA_CHECK(cudaFree(d_x));
    CUDA_CHECK(cudaFree(d_q));
    CUDA_CHECK(cudaFree(d_scalar));
  }  



  // void act_gpu2( T* d_res, const T* d_v ) const {
  //   assert( on_gpu );

  //   cusparseHandle_t handle = NULL;

  //   const cusparseDirection_t dirA = CUSPARSE_DIRECTION_COLUMN;
  //   constexpr Idx blockDim = ComilationConst::N;
  //   constexpr int mb = 1;
  //   constexpr int nb = 1;
  //   constexpr int nnzb = 1;

  //   int* bsrRowPtr;
  //   cudaMalloc((void**)&bsrRowPtrC, sizeof(int) *(mb+1));

  //   CuC* descrA
  //   cusparseXcsr2bsrNnz(handle, dirA, m, n,

  // 			descrA, csrRowPtrA, csrColIndA, blockDim,
  // 			descrC, bsrRowPtrC, &nnzb);

  //   cudaMalloc((void**)&bsrColIndC, sizeof(int)*nnzb);
  //   cudaMalloc((void**)&bsrValC, sizeof(float)*(blockDim*blockDim)*nnzb);
  //   cusparseScsr2bsr(handle, dirA, m, n,
  // 		     descrA, csrValA, csrRowPtrA, csrColIndA, blockDim,
  // 		     descrC, bsrValC, bsrRowPtrC, bsrColIndC);
  //   // step 2: allocate vector x and vector y large enough for bsrmv
  //   cudaMalloc((void**)&x, sizeof(float)*(nb*blockDim));
  //   cudaMalloc((void**)&y, sizeof(float)*(mb*blockDim));
  //   cudaMemcpy(x, hx, sizeof(float)*n, cudaMemcpyHostToDevice);
  //   cudaMemcpy(y, hy, sizeof(float)*m, cudaMemcpyHostToDevice);
  //   // step 3: perform bsrmv
  //   cusparseSbsrmv(handle, dirA, transA, mb, nb, nnzb, &alpha,
  // 		   descrC, bsrValC, bsrRowPtrC, bsrColIndC, blockDim, x, &beta, y);

  //     cusparseZbsrmv(cusparseHandle_t         handle,
  //              cusparseDirection_t      dir,
  //              cusparseOperation_t      trans,
  //              int                      mb,
  //              int                      nb,
  //              int                      nnzb,
  //              const cuDoubleComplex*   alpha,
  //              const cusparseMatDescr_t descr,
  //              const cuDoubleComplex*   bsrVal,
  //              const int*               bsrRowPtr,
  //              const int*               bsrColInd,
  //              int                      blockDim,
  //              const cuDoubleComplex*   x,
  //              const cuDoubleComplex*   beta,
  //              cuDoubleComplex*         y)


  //   constexpr Idx N = CompilationConst::N;
  //   mult<T,N><<<NBlocks, NThreadsPerBlock>>>(d_res,
  // 					     d_v,
  // 					     this->val,
  // 					     this->cols,
  // 					     this->rows);
  // }





  using Link = std::array<Idx,2>; // <int,int>;

  const Idx ix=2;
  const Idx iy=lattice.nns[ix][0];
  const Link ell{ix,iy};

  std::vector<COOEntry> en;
  DW.d_coo_format(en, U, ell);

  std::vector<Complex> Dxi(N), xi(N);
  for(int i=0; i<N; i++) {
    xi[i] = rng.gaussian();
  }
  matmulcoo<N>(reinterpret_cast<CuC*>(Dxi.data()),
               reinterpret_cast<CuC*>(xi.data()),
               en
               );

  std::cout << "Dxi = " << std::endl;
  for(auto elem:Dxi) std::cout << elem << std::endl;
  // std::cout << "is = " << std::endl;
  // for(auto elem:is) std::cout << elem << std::endl;
  // std::cout << "js = " << std::endl;
  // for(auto elem:js) std::cout << elem << std::endl;
  {
    // for(auto elem:en) std::cout << elem << std::endl;
    std::sort( en.begin(), en.end() );
    // for(auto elem:en) std::cout << elem << std::endl;

    Idx len=en.size();
    std::vector<CuC> v(len);
    std::vector<Idx> cols(len);
    std::vector<Idx> rows;

    for(Idx k=0; k<len; k++){
      v[k] = en[k].v;
      cols[k] = en[k].j;
    }

    Idx k=0;
    rows.push_back(k);
    for(Idx i=0; i<N; i++){
      while(en[k].i == i) k++;
      rows.push_back(k);
    }

    matmul<N>(reinterpret_cast<CuC*>(Dxi.data()),
              reinterpret_cast<CuC*>(xi.data()),
              v, cols, rows
              );

    std::cout << "Dxi = " << std::endl;
    for(auto elem:Dxi) std::cout << elem << std::endl;


    // rows.push_back(N);

    // std::cout << "v = " << std::endl;
    // for(auto elem:v) std::cout << real(elem) << " " << imag(elem) << std::endl;
    // std::cout << "is = " << std::endl;
    // for(auto elem:cols) std::cout << elem << std::endl;
    // std::cout << "js = " << std::endl;
    // for(auto elem:rows) std::cout << elem << std::endl;

    // std::cout << v.size() << " "
    //           << cols.size() << " "
    //           << rows.size() << std::endl;

    // std::cout << len << " " << N << std::endl;
  }


  // bool on_gpu = false;
  // T* val;
  // Idx* cols;
  // Idx* rows;
  // SparseMatrix();

  // void act_cpu( std::vector<T>& res, const std::vector<T>& v ) const {
  //   assert( !on_gpu );

  //   constexpr Idx N = CompilationConst::N;
  //   for(Idx i=0; i<N; i++) {
  //     res[i] = cplx(0.0);
  //     const int row_start = rows[i];
  //     const int row_end = rows[i+1];
  //     for(int jj=row_start; jj<row_end; jj++) res[i] = res[i] + val[jj] * v[ cols[jj] ];
  //   }
  // }


  COO coo;
  DW.d_coo_format(coo.en, U, ell);
  coo.set();

  std::vector<Complex> Dxi(N), xi(N);
  for(int i=0; i<N; i++) {
    xi[i] = rng.gaussian();
  }

  {
    CuC *d_Dxi, *d_xi;
    CUDA_CHECK(cudaMalloc(&d_Dxi, N*CD));
    CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
    CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

    coo( d_Dxi, d_xi );

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Dxi.data()), d_Dxi, N*CD, D2H));

    CUDA_CHECK(cudaFree(d_Dxi));
    CUDA_CHECK(cudaFree(d_xi));
  }

  std::cout << "Dxi = " << std::endl;
  for(auto elem:Dxi) std::cout << elem << std::endl;



  //   void d_coo_format( std::vector<Complex>& v,
  //       	     std::vector<Idx>& is,
  //       	     std::vector<Idx>& js,
  //       	     const Gauge& U,
  //       	     const Link& ell ) const {
  //   const Idx ix0 = ell[0];
  //   const Idx iy0 = ell[1];

  //   v.clear();
  //   is.clear();
  //   js.clear();

  //   {
  //     // pos
  //     const Idx ix = ix0;
  //     for(int jj=0; jj<3; jj++){
  //       const Idx iy = lattice.nns[ix][jj];
  //       if(iy!=iy0) continue;
  //       const MS tmp = lattice.vol[ix]/lattice.mean_vol * (r*lattice.u[ix][jj]*sigma[0] - gamma(ix, jj)) * I*std::exp( I* U(Link{ix,iy})) * Omega(ix, iy);

  //       // res[NS*ix] += -tmp(0,0)*v[NS*iy] - tmp(0,1)*v[NS*iy+1];
  //       v.push_back(-tmp(0,0)); is.push_back(NS*ix); js.push_back(NS*iy);
  //       v.push_back(-tmp(0,1)); is.push_back(NS*ix); js.push_back(NS*iy+1);

  //       // res[NS*ix+1] += -tmp(1,0)*v[NS*iy] - tmp(1,1)*v[NS*iy+1];
  //       v.push_back(-tmp(1,0)); is.push_back(NS*ix+1); js.push_back(NS*iy);
  //       v.push_back(-tmp(1,1)); is.push_back(NS*ix+1); js.push_back(NS*iy+1);
  //     }
  //   }

  //   {
  //     // neg
  //     const Idx iy = iy0;
  //     for(int jj=0; jj<3; jj++){
  //       const Idx ix = lattice.nns[iy0][jj];
  //       if(ix!=ix0) continue;
  //       const MS tmp = -lattice.vol[iy]/lattice.mean_vol * (r*lattice.u[iy][jj]*sigma[0] - gamma(iy, jj)) * I*std::exp( I* U(Link{iy,ix})) * Omega(iy, ix);

  //       // res[NS*iy] += -tmp(0,0)*v[NS*ix] - tmp(0,1)*v[NS*ix+1];
  //       v.push_back(-tmp(0,0)); is.push_back(NS*iy); js.push_back(NS*ix);
  //       v.push_back(-tmp(0,1)); is.push_back(NS*iy); js.push_back(NS*ix+1);

  //       // res[NS*iy+1] += -tmp(1,0)*v[NS*ix] - tmp(1,1)*v[NS*ix+1];
  //       v.push_back(-tmp(1,0)); is.push_back(NS*iy+1); js.push_back(NS*ix);
  //       v.push_back(-tmp(1,1)); is.push_back(NS*iy+1); js.push_back(NS*ix+1);
  //     }
  //   }
  // }


  {
    Overlap Dov(DW);
    Dov.compute(U);

    std::vector<Complex> xi(N), Dxi(N), DHDxi(N), DHDxi2(N);
    for(int i=0; i<N; i++) xi[i] = rng.gaussian();

    Dov( Dxi, xi );
    Dov.adj( DHDxi, Dxi );
    Dov.sq( DHDxi2, xi );

    for(int i=0; i<N; i++) {
      std::cout << xi[i] << " " << DHDxi[i] << " " << DHDxi2[i] << std::endl;
    }

  }



  // {
  //   Overlap Dov(DW);
  //   Dov.compute(U);

  //   std::vector<Complex> xi(N), Dxi(N), DHDxi(N), DHDxi2(N);
  //   for(int i=0; i<N; i++) xi[i] = rng.gaussian();

  //   Dov.mult( Dxi, xi );
  //   Dov.adj( DHDxi, Dxi );
  //   Dov.sq( DHDxi2, xi );

  //   for(int i=0; i<N; i++) {
  //     std::cout << xi[i] << " " << DHDxi[i] << " " << DHDxi2[i] << std::endl;
  //   }
  // }


  // {
  //   Overlap Dov(DW);
  //   Dov.compute(U);

  //   std::vector<Complex> xi(N), Dxi(N), Dxi2(N);
  //   for(int i=0; i<N; i++) xi[i] = rng.gaussian();

  //   Dov.mult( Dxi, xi );

  //   CuC *d_Dxi2, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Dxi2, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   Dov.mult( d_Dxi2, d_xi );

  //   CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Dxi2.data()), d_Dxi2, N*CD, D2H));

  //   CUDA_CHECK(cudaFree(d_Dxi2));
  //   CUDA_CHECK(cudaFree(d_xi));


  //   for(int i=0; i<N; i++) {
  //     std::cout << xi[i] << " " << Dxi[i] << " " << Dxi2[i] << std::endl;
  //   }
  // }


  // {
  //   Overlap Dov(DW);
  //   Dov.compute(U);

  //   std::vector<Complex> xi(N), Dxi(N), Dxi2(N);
  //   for(int i=0; i<N; i++) xi[i] = rng.gaussian();

  //   Dov.adj( Dxi, xi );

  //   CuC *d_Dxi2, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Dxi2, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   Dov.adj( d_Dxi2, d_xi );

  //   CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Dxi2.data()), d_Dxi2, N*CD, D2H));

  //   CUDA_CHECK(cudaFree(d_Dxi2));
  //   CUDA_CHECK(cudaFree(d_xi));


  //   for(int i=0; i<N; i++) {
  //     std::cout << xi[i] << " " << Dxi[i] << " " << Dxi2[i] << std::endl;
  //   }
  // }


  // {
  //   Overlap Dov(DW);
  //   Dov.compute(U);

  //   std::vector<Complex> xi(N), Dxi(N), Dxi2(N);
  //   for(int i=0; i<N; i++) xi[i] = rng.gaussian();

  //   Dov.sq( Dxi, xi );

  //   CuC *d_Dxi2, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Dxi2, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   Dov.sq( d_Dxi2, d_xi );

  //   CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Dxi2.data()), d_Dxi2, N*CD, D2H));

  //   CUDA_CHECK(cudaFree(d_Dxi2));
  //   CUDA_CHECK(cudaFree(d_xi));


  //   for(int i=0; i<N; i++) {
  //     std::cout << xi[i] << " " << Dxi[i] << " " << Dxi2[i] << std::endl;
  //   }
  // }



  {
    Overlap Dov(DW);
    Dov.compute(U);

    std::vector<Complex> xi(N), Dxi(N), Dxi2(N);
    for(int i=0; i<N; i++) xi[i] = rng.gaussian();

    Dov.mult( Dxi, xi );

    CuC *d_Dxi2, *d_xi;
    CUDA_CHECK(cudaMalloc(&d_Dxi2, N*CD));
    CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
    CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

    // Dov.mult_device( d_Dxi2, d_xi );
    auto f = std::bind(&Overlap::mult_device, &Dov, std::placeholders::_1, std::placeholders::_2);
    // f( d_Dxi2, d_xi );
    // std::cout << typeid(f).name() << std::endl;

    // auto f = std::bind<LinOpWrapper::Function**, LinOpWrapper, CuC*, CuC*>(&Overlap::mult, Dov, std::placeholders::_1, std::placeholders::_2);
    // LinOpWrapper Op( std::bind(&Overlap::mult_device, &Dov, std::placeholders::_1, std::placeholders::_2) );
    LinOpWrapper Op( f );
    Op( d_Dxi2, d_xi );

    std::vector<Complex> eta(N);
    {
      MatPoly Op;
      Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );
      Op.solve<N>( eta, xi );
    }


    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(Dxi2.data()), d_Dxi2, N*CD, D2H));
    CUDA_CHECK(cudaFree(d_Dxi2));
    CUDA_CHECK(cudaFree(d_xi));

    for(int i=0; i<N; i++) {
      std::cout << xi[i] << " " << Dxi[i] << " " << Dxi2[i] << std::endl;
    }
  }



  {
    Overlap Dov(DW);
    Dov.compute(U);

    std::vector<Complex> xi(N), xi2(N); // , Dxi(N), Dxi2(N);
    for(int i=0; i<N; i++) xi[i] = rng.gaussian();

    // Dov.mult( Dxi, xi );

    CuC *d_xi, *d_eta;
    CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
    CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
    CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

    auto f = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
    LinOpWrapper Op( f );

    // for(int i=0; i<100; i++) Op( d_eta, d_xi );

    std::vector<Complex> eta(N);
    {
      MatPoly Poly;
      Poly.push_back ( cplx(1.0), {&Op} );
      Poly.solve<N>( d_eta, d_xi );
    }

    f( d_xi, d_eta );

    // CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(eta.data()), d_eta, N*CD, D2H));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xi2.data()), d_xi, N*CD, D2H));
    CUDA_CHECK(cudaFree(d_xi));
    CUDA_CHECK(cudaFree(d_eta));

    for(int i=0; i<N; i++) {
      std::cout << xi[i] - xi2[i] << std::endl;
    }
  }




  // {
  //   Dov.compute(U);

  //   CuC *d_Z, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Z, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   int m=1;
  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0/(Dov.lambda_max*Dov.lambda_max)), {&Dov.M_DW, &Dov.M_DWH} );
  //   const CuC a = cplx(-Dov.k*Dov.k/Dov.cp[m]);
  //   Op.push_back ( a, {} );
  //   Op.solve<N>( d_Z, d_xi );
  //   Op.Zdscal<N>( Dov.A[m], d_Z );

  //   Dummy.dot<N>( &Sf, d_xi, d_Z );
  //   std::cout << "S = " << real(Sf) << " " << imag(Sf) << std::endl;

  //   CUDA_CHECK(cudaFree(d_Z));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }

  // {
  //   Dov.compute(U);

  //   CuC *d_Z, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Z, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   int m=1;
  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0/(Dov.lambda_max*Dov.lambda_max)), {&Dov.M_DW, &Dov.M_DWH} );
  //   const CuC a = cplx(-Dov.k*Dov.k/Dov.cp[m]);
  //   Op.push_back ( a, {} );
  //   Op.solve<N>( d_Z, d_xi );

  //   // --------------

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, ell);
  //   coo.do_it();

  //   CuC *d_tmp1, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

  //   coo( d_tmp1, d_Z ); // DH
  //   Dov.M_DWH( d_tmp2, d_tmp1 );

  //   Op.Zdscal<N>( Dov.A[m], d_Z );

  //   Dummy.dot<N>( &grad, d_Z, d_tmp2 );

  //   CUDA_CHECK(cudaFree(d_tmp1));
  //   CUDA_CHECK(cudaFree(d_tmp2));

  //   // --------------

  //   grad = -2.0 * grad * cplx(1.0/(Dov.lambda_max*Dov.lambda_max));
  //   std::cout << "grad = " << real(grad) << " " << imag(grad) << std::endl;

  //   CUDA_CHECK(cudaFree(d_Z));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }


  // {
  //   Dov.compute(UP);

  //   CuC *d_Z, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Z, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   int m=1;
  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0/(Dov.lambda_max*Dov.lambda_max)), {&Dov.M_DW, &Dov.M_DWH} );
  //   const CuC a = cplx(-Dov.k*Dov.k/Dov.cp[m]);
  //   Op.push_back ( a, {} );
  //   Op.solve<N>( d_Z, d_xi );
  //   Op.Zdscal<N>( Dov.A[m], d_Z );

  //   Dummy.dot<N>( &Sfp, d_xi, d_Z );
  //   std::cout << "Sp = " << real(Sfp) << " " << imag(Sfp) << std::endl;

  //   CUDA_CHECK(cudaFree(d_Z));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }


  // {
  //   Dov.compute(UM);

  //   CuC *d_Z, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Z, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   int m=1;
  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0/(Dov.lambda_max*Dov.lambda_max)), {&Dov.M_DW, &Dov.M_DWH} );
  //   const CuC a = cplx(-Dov.k*Dov.k/Dov.cp[m]);
  //   Op.push_back ( a, {} );
  //   Op.solve<N>( d_Z, d_xi );
  //   Op.Zdscal<N>( Dov.A[m], d_Z );

  //   Dummy.dot<N>( &Sfm, d_xi, d_Z );
  //   std::cout << "Sm = " << real(Sfm) << " " << imag(Sfm) << std::endl;

  //   CUDA_CHECK(cudaFree(d_Z));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }

  // CuC check = (Sfp-Sfm)/(2.0*eps);
  // std::cout << "grad = " << real(grad) << " " << imag(grad) << std::endl;
  // std::cout << "check = " << real(check) << " " << imag(check) << std::endl;



  // std::cout << std::endl << std::endl << std::endl << std::endl;




  // // auto f_Dov = std::bind(&Overlap::mult_device2, &Dov, std::placeholders::_1, std::placeholders::_2);
  // // auto f_Dov = std::bind(&Overlap::mult_device3, &Dov, std::placeholders::_1, std::placeholders::_2);
  // auto f_Dov = std::bind(&Overlap::mult_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  // LinOpWrapper M_Dov( f_Dov );

  // {
  //   Dov.compute(U);

  //   CuC *d_eta, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   M_Dov( d_eta, d_xi );

  //   Dummy.dot<N>( &Sf, d_xi, d_eta );
  //   std::cout << "S = " << real(Sf) << " " << imag(Sf) << std::endl;

  //   CUDA_CHECK(cudaFree(d_eta));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }

  // {
  //   Dov.compute(U);

  //   CuC *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));
  //   // grad_d = Dov.grad2( ell, U, d_xi ); // DH
  //   grad_d = Dov.grad3( ell, U, d_xi ); // DH
  //   CUDA_CHECK(cudaFree(d_xi));
  // }

  // {
  //   Dov.compute(UP);

  //   CuC *d_eta, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   M_Dov( d_eta, d_xi );

  //   Dummy.dot<N>( &Sfp, d_xi, d_eta );
  //   std::cout << "Sp = " << real(Sfp) << " " << imag(Sfp) << std::endl;

  //   CUDA_CHECK(cudaFree(d_eta));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }

  // {
  //   Dov.compute(UM);

  //   CuC *d_eta, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   M_Dov( d_eta, d_xi );

  //   Dummy.dot<N>( &Sfm, d_xi, d_eta );
  //   std::cout << "Sm = " << real(Sfm) << " " << imag(Sfm) << std::endl;

  //   CUDA_CHECK(cudaFree(d_eta));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }


  // CuC check2 = (Sfp-Sfm)/(2.0*eps);
  // std::cout << "grad = " << 0.5 * grad_d << std::endl;
  // // std::cout << "grad = " << real(grad) << " " << imag(grad) << std::endl;
  // std::cout << "check = " << real(check2) << " " << imag(check2) << std::endl;







  // {
  //   CuC *d_Dxi, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_Dxi, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   COO coo;
  //   MatPoly Dummy;
  //   DW.d_coo_format(coo.en, U, ell);
  //   CuC dS;

  //   coo.do_it();
  //   coo( d_Dxi, d_xi );
  //   Dummy.dot<N>( &dS, d_xi, d_Dxi );
  //   std::cout << "dS = " << real(dS) << " " << imag(dS) << std::endl;

  //   coo.do_conjugate();
  //   coo( d_Dxi, d_xi ); // DH
  //   Dummy.dot<N>( &dS, d_Dxi, d_xi );
  //   std::cout << "dS = " << real(dS) << " " << imag(dS) << std::endl;

  //   // -----------------------

  //   CSR M_DW;
  //   DWDevice d_DW(DW); // actual data used in M_DW, M_DWH
  //   d_DW.associateCSR( M_DW, false );

  //   CuC Sfp, Sfm;
  //   const double eps = 1.0e-5;

  //   Gauge UP(U);
  //   UP[il] += eps;
  //   d_DW.update( UP );
  //   M_DW( d_Dxi, d_xi );
  //   Dummy.dot<N>( &Sfp, d_xi, d_Dxi );

  //   Gauge UM(U);
  //   UM[il] -= eps;
  //   d_DW.update( UM );
  //   M_DW( d_Dxi, d_xi );
  //   Dummy.dot<N>( &Sfm, d_xi, d_Dxi );

  //   CUDA_CHECK(cudaFree(d_Dxi));
  //   CUDA_CHECK(cudaFree(d_xi));

  //   CuC num = (Sfp-Sfm)/(2.0*eps);
  //   std::cout << "dS = " << real(num) << " " << imag(num) << std::endl;
  // }



  // const Idx ix=2;
  // const Idx iy=lattice.nns[ix][0];
  // const Link ell{ix,iy};

  // COO coo;
  // DW.d_coo_format(coo.en, U, ell);
  // coo.do_it();




  // {
  //   Overlap Dov(DW);
  //   Dov.compute(U);

  //   std::vector<Complex> xi(N), xi2(N); // , Dxi(N), Dxi2(N);
  //   for(int i=0; i<N; i++) xi[i] = rng.gaussian();

  //   // Dov.mult( Dxi, xi );

  //   CuC *d_xi, *d_eta;
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   auto f = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  //   LinOpWrapper Op( f );

  //   // for(int i=0; i<100; i++) Op( d_eta, d_xi );

  //   std::vector<Complex> eta(N);
  //   {
  //     MatPoly Poly;
  //     Poly.push_back ( cplx(1.0), {&Op} );
  //     Poly.solve<N>( d_eta, d_xi );
  //   }

  //   f( d_xi, d_eta );

  //   // CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(eta.data()), d_eta, N*CD, D2H));
  //   CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xi2.data()), d_xi, N*CD, D2H));
  //   CUDA_CHECK(cudaFree(d_xi));
  //   CUDA_CHECK(cudaFree(d_eta));

  //   for(int i=0; i<N; i++) {
  //     std::cout << xi[i] << " " << xi2[i] << " " << xi[i] - xi2[i] << std::endl;
  //   }
  // }





  //   auto f_Dov = std::bind(&Overlap::mult_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  // LinOpWrapper M_Dov( f_Dov );


  // {
  //   Dov.compute(U);

  //   CuC *d_eta, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0), {&M_Dov} );
  //   Op.solve<N>( d_eta, d_xi );

  //   Dummy.dot<N>( &Sf, d_xi, d_eta );
  //   std::cout << "S = " << real(Sf) << " " << imag(Sf) << std::endl;

  //   CUDA_CHECK(cudaFree(d_eta));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }

  // {
  //   Dov.compute(U);

  //   CuC *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));
  //   grad = Dov.grad( ell, U, d_xi ); // DH

  //   // --------------

  //   // std::cout << "grad = " << real(grad) << " " << imag(grad) << std::endl;
  //   std::cout << "grad = " << grad << std::endl;

  //   CUDA_CHECK(cudaFree(d_xi));
  // }



  // {
  //   Dov.compute(UP);

  //   CuC *d_eta, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0), {&M_Dov} );
  //   Op.solve<N>( d_eta, d_xi );

  //   Dummy.dot<N>( &Sfp, d_xi, d_eta );
  //   std::cout << "Sp = " << real(Sfp) << " " << imag(Sfp) << std::endl;

  //   CUDA_CHECK(cudaFree(d_eta));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }

  // {
  //   Dov.compute(UM);

  //   CuC *d_eta, *d_xi;
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0), {&M_Dov} );
  //   Op.solve<N>( d_eta, d_xi );

  //   Dummy.dot<N>( &Sfm, d_xi, d_eta );
  //   std::cout << "Sm = " << real(Sfm) << " " << imag(Sfm) << std::endl;

  //   CUDA_CHECK(cudaFree(d_eta));
  //   CUDA_CHECK(cudaFree(d_xi));
  // }


  // CuC check = (Sfp-Sfm)/(2.0*eps);
  // std::cout << "grad = " << grad << std::endl;
  // // std::cout << "grad = " << real(grad) << " " << imag(grad) << std::endl;
  // std::cout << "check = " << real(check) << " " << imag(check) << std::endl;






  




  //   COO coo;
  //   MatPoly Dummy;
  //   DW.d_coo_format(coo.en, U, ell);
  //   CuC dS;


  //   coo.do_it();
  //   coo( d_Dxi, d_xi );
  //   Dummy.dot<N>( &dS, d_xi, d_Dxi );
  //   std::cout << "dS = " << real(dS) << " " << imag(dS) << std::endl;

  //   coo.do_conjugate();
  //   coo( d_Dxi, d_xi ); // DH
  //   Dummy.dot<N>( &dS, d_Dxi, d_xi );
  //   std::cout << "dS = " << real(dS) << " " << imag(dS) << std::endl;

  //   // -----------------------

  //   CSR M_DW;
  //   DWDevice d_DW(DW); // actual data used in M_DW, M_DWH
  //   d_DW.associateCSR( M_DW, false );

  //   CuC Sfp, Sfm;

  //   d_DW.update( UP );
  //   M_DW( d_Dxi, d_xi );
  //   Dummy.dot<N>( &Sfp, d_xi, d_Dxi );

  //   d_DW.update( UM );
  //   M_DW( d_Dxi, d_xi );
  //   Dummy.dot<N>( &Sfm, d_xi, d_Dxi );

  //   CUDA_CHECK(cudaFree(d_Dxi));
  //   CUDA_CHECK(cudaFree(d_xi));

  //   CuC num = (Sfp-Sfm)/(2.0*eps);
  //   std::cout << "dS = " << real(num) << " " << imag(num) << std::endl;
  // }
















  // {
  //   Overlap Dov(DW, 0.0001, 21);
  //   Dov.compute(U);

  //   // std::vector<Complex> xi(N), eta(N);
  //   // for(int i=0; i<N; i++) xi[i] = rng.gaussian();

  //   // -------------

  //   CuC *d_xi, *d_eta, *d_dummy;
  //   CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_dummy, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  //   // auto f_DovH_Dov = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  //   auto f_DovH_Dov = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  //   LinOpWrapper M_DovH_Dov( f_DovH_Dov );

  //   double TOL=1.0e-12;
  //   // double Sf = 0.0;
  //   // {
  //   //   MatPoly Op;
  //   //   Op.push_back ( cplx(1.0), {&M_DovH_Dov} );

  //   //   // Op.solve<N>( d_eta, d_xi );
  //   //   Op.on_gpu<N>( d_eta, d_xi );

  //   //   Complex tmp;
  //   //   Op.dot<N>( reinterpret_cast<CuC*>(&tmp), d_xi, d_eta );
  //   //   assert( std::abs(tmp.imag()/tmp.real()) < TOL );
  //   //   Sf = tmp.real();
  //   // }

  //   // // CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(xi.data()), d_xi, N*CD, D2H));
  //   // // CUDA_CHECK(cudaMemcpy(reinterpret_cast<CuC*>(eta.data()), d_eta, N*CD, D2H));
  //   // // for(int i=0; i<N; i++) {free
  //   // //   std::cout << xi[i] << " " << eta[i] << std::endl;
  //   // // }
  //   // std::cout << Sf << std::endl;

  //   // -------------

  //   // Idx il=0;

  //   std::cout << "######## calculating grad. " << std::endl;
  //   double dSf = Dov.grad3( lattice.links[il], U, d_xi );
  //   // std::cout << "grad = " << dSf << std::endl;

  //   // std::cout << "######## calculating mult. " << std::endl;
  //   // Dov.mult_device( d_dummy, d_phi );
  //   // std::cout << "######## . " << std::endl;

  //   // -------------

  //   double Sfp = 0.0, Sfm = 0.0;
  //   const double eps = 1.0e-5;

  //   {
  //     {
  //       // Gauge UP(U);
  //       // UP[il] += eps;
  //       Dov.compute(UP);

  //       MatPoly Op;
  //       // Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );
  //       // Op.solve<N>( eta, xi );
  //       // for(Idx i=0; i<N; i++) Sfp += std::real( std::conj(xi[i]) * eta[i] );
  //       Op.push_back ( cplx(1.0), {&M_DovH_Dov} );
  //       // Op.solve<N>( d_eta, d_xi );
  //       Op.on_gpu<N>( d_eta, d_xi );

  //       Complex tmp;
  //       Op.dot<N>( reinterpret_cast<CuC*>(&tmp), d_xi, d_eta );
  //       assert( std::abs(tmp.imag()/tmp.real()) < TOL );
  //       Sfp = tmp.real();
  //     }

  //     {
  //       Gauge UM(U);
  //       UM[il] -= eps;
  //       Dov.compute(UM);

  //       MatPoly Op;
  //       // Op.push_back ( cplx(1.0), {&(Dov.M_DW), &(Dov.M_DWH)} );
  //       // Op.solve<N>( eta, xi );
  //       // for(Idx i=0; i<N; i++) Sfm += std::real( std::conj(xi[i]) * eta[i] );
  //       Op.push_back ( cplx(1.0), {&M_DovH_Dov} );
  //       // Op.solve<N>( d_eta, d_xi );
  //       Op.on_gpu<N>( d_eta, d_xi );

  //       Complex tmp;
  //       Op.dot<N>( reinterpret_cast<CuC*>(&tmp), d_xi, d_eta );
  //       assert( std::abs(tmp.imag()/tmp.real()) < TOL );
  //       Sfm = tmp.real();
  //     }

  //   }

  //   std::cout << "grad = " << dSf << std::endl;
  //   std::cout << "check = " << (Sfp-Sfm)/(2.0*eps) << std::endl;

  //   CUDA_CHECK(cudaFree(d_xi));
  //   CUDA_CHECK(cudaFree(d_eta));
  //   CUDA_CHECK(cudaFree(d_dummy));
  // }



  Idx il=0;
  Link ell = lattice.links[il];

  std::vector<Complex> Dxi(N), xi(N);
  for(int i=0; i<N; i++) xi[i] = rng.gaussian() + 1.0*Complex(0.0,1.0)*rng.gaussian();


  Overlap Dov(DW, 0.0001, 11);

  // CuC Sf, Sfp, Sfm;
  CuC Sf, Sfp, Sfm, grad;
  double grad_d;

  MatPoly Dummy;
  CuC dS;

  const double eps = 1.0e-5;
  Gauge UP(U);
  UP[il] += eps;
  Gauge UM(U);
  UM[il] -= eps;

  // auto f_Dov = std::bind(&Overlap::mult_device2, &Dov, std::placeholders::_1, std::placeholders::_2);
  // auto f_Dov = std::bind(&Overlap::mult_device3, &Dov, std::placeholders::_1, std::placeholders::_2);
  // auto f_Dov = std::bind(&Overlap::mult_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  auto f_DHDov = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DHDov( f_DHDov );

  MatPoly OpDHDov;
  OpDHDov.push_back ( cplx(1.0), {&M_DHDov} );

  CuC *d_eta, *d_xi;
  CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  {
    Dov.compute(U);
    OpDHDov.solve<N>( d_eta, d_xi );
    OpDHDov.dot<N>( &Sf, d_xi, d_eta );
    std::cout << "S = " << real(Sf) << " " << imag(Sf) << std::endl;
  }

  {
    Dov.compute(U);
    grad_d = Dov.grad_device( ell, U, d_eta ); // DH
  }

  {
    Dov.compute(UP);
    OpDHDov.solve<N>( d_eta, d_xi );
    OpDHDov.dot<N>( &Sfp, d_xi, d_eta );
    std::cout << "Sp = " << real(Sfp) << " " << imag(Sfp) << std::endl;
  }

  {
    Dov.compute(UM);
    OpDHDov.solve<N>( d_eta, d_xi );
    OpDHDov.dot<N>( &Sfm, d_xi, d_eta );
    std::cout << "Sm = " << real(Sfm) << " " << imag(Sfm) << std::endl;
  }

  CuC check2 = (Sfp-Sfm)/(2.0*eps);
  std::cout << "grad = " << grad_d << std::endl;
  // std::cout << "grad = " << real(grad) << " " << imag(grad) << std::endl;
  std::cout << "check = " << real(check2) << " " << imag(check2) << std::endl;


  CUDA_CHECK(cudaFree(d_eta));
  CUDA_CHECK(cudaFree(d_xi));








  // void adj(std::vector<Complex>& res, const std::vector<Complex>& xi) const {
  //   std::vector<Complex> DHxi(xi.size());

  //   {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //     Op.from_cpu<N>( DHxi, xi );
  //     for(Idx i=0; i<res.size(); i++) DHxi[i] *= C;
  //   }

  //   res = DHxi;

  //   std::vector<Complex> tmp(xi.size());
  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( tmp, DHxi );
  //     for(Idx i=0; i<res.size(); i++) res[i] += A[m] * tmp[i];
  //   }

  //   for(Idx i=0; i<res.size(); i++) res[i] += xi[i];
  // }


  // void sq(std::vector<Complex>& res, const std::vector<Complex>& xi) const {
  //   std::vector<Complex> tmp1(xi.size()), tmp2(xi.size());
  //   this->mult(tmp1, xi);
  //   this->adj(tmp2, xi);
  //   for(Idx i=0;  i<res.size(); i++) res[i] = tmp1[i] + tmp2[i];
  // }



  // double grad( const Link& link, const Gauge& U, const std::vector<Complex>& eta ) const {
  //   double res = 0.0;

  //   std::vector<Complex> Xdag_eta(eta.size());
  //   {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //     Op.from_cpu<N>(Xdag_eta, eta);
  //   }

  //   std::vector<std::vector<Complex>> Zs(size, std::vector<Complex>(eta.size()) );
  //   std::vector<std::vector<Complex>> Ys(size, std::vector<Complex>(eta.size()) );

  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( Zs[m], eta );
  //     Op.solve<N>( Ys[m], Xdag_eta );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   {
  //     std::vector<Complex> sum(eta.size(), 0.0);
  //     for(int m=1; m<size; m++) {
  //       for(Idx i=0; i<eta.size(); i++) sum[i] += Zs[m][i];
  //     }
  //     std::vector<Complex> dD_sum(eta.size(), 0.0);
  //     matmulcoo<N>( reinterpret_cast<CuC*>(dD_sum.data()),
  //                   reinterpret_cast<const CuC*>(sum.data()),
  //                    coo.en );
  //     for(Idx i=0; i<eta.size(); i++) res += std::real( std::conj(eta[i]) * dD_sum[i] );
  //   }
  //   {
  //     for(int m=1; m<size; m++) {
  //       std::vector<Complex> XZm(eta.size(), 0.0);
  //       {
  //         MatPoly Op;
  //         Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //         Op.from_cpu<N>(XZm, Zs[m]);
  //       }
  //       std::vector<Complex> dD_Ym(eta.size(), 0.0);
  //       matmulcoo<N>( reinterpret_cast<CuC*>(dD_Ym.data()),
  //                     reinterpret_cast<const CuC*>(Ys[m].data()),
  //                     coo.en );
  //       for(Idx i=0; i<eta.size(); i++) res -= std::real( std::conj(XZm[i]) * dD_Ym[i] );
  //     }
  //   }
  //   {
  //     for(int m=1; m<size; m++) {
  //       std::vector<Complex> XYm(eta.size(), 0.0);
  //       {
  //         MatPoly Op;
  //         Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //         Op.from_cpu<N>(XYm, Ys[m]);
  //       }
  //       std::vector<Complex> dD_Zm(eta.size(), 0.0);
  //       matmulcoo<N>( reinterpret_cast<CuC*>(dD_Zm.data()),
  //                     reinterpret_cast<const CuC*>(Zs[m].data()),
  //                     coo.en );
  //       for(Idx i=0; i<eta.size(); i++) res -= std::real( std::conj(XYm[i]) * dD_Zm[i] );
  //     }
  //   }

  //   res *= -2.0 * E * 2.0 / (1.0+lambda_inv) / (k*M);

  //   return res;
  // }


  // void mult(std::vector<Complex>& res, const std::vector<Complex>& xi) const {
  //   res = xi;
  //   std::vector<Complex> tmp(xi.size());

  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );
  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( tmp, xi );
  //     for(Idx i=0; i<res.size(); i++) res[i] += A[m] * tmp[i];
  //   }

  //   for(Idx i=0; i<res.size(); i++) tmp[i] = E * 2.0 / (1.0+lambda_inv) / (k*M) * res[i];
  //   MatPoly Op;
  //   Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   Op.from_cpu<N>( res, tmp );

  //   for(Idx i=0; i<res.size(); i++) res[i] += xi[i];
  // }

  // MatPoly Dummy;
  // double mu;
  // Dummy.dot2self<N>(&mu, d_xi);
  // std::cout << "|d_xi| = " << mu << std::endl;



  // double grad( const Link& link, const Gauge& U, CuC* d_eta ) const {
  //   double res = 0.0;

  //   double mu;
  //   CuC inner;
  //   MatPoly Dummy;

  //   CuC* d_Xdag_eta;
  //   CUDA_CHECK(cudaMalloc(&d_Xdag_eta, N*CD));
  //   CUDA_CHECK(cudaMemcpy(d_Xdag_eta, d_eta, N*CD, D2D)); // @@@@@@
  //   // {
  //   //   MatPoly Op;
  //   //   Op.push_back ( cplx(1.0/(lambda_max)), {&M_DWH} );
  //   //   Op.on_gpu<N>(d_Xdag_eta, d_eta);
  //   // }

  //   std::vector<CuC*> d_Zs(size);
  //   std::vector<CuC*> d_Ys(size);
  //   for(int m=1; m<size; m++) {
  //     CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));
  //     CUDA_CHECK(cudaMalloc(&d_Ys[m], N*CD));
  //   }

  //   for(int m=1; m<size; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Zs[m], d_eta );
  //     Op.solve<N>( d_Ys[m], d_Xdag_eta );

  //     Op.Zdscal<N>( A[m], d_Zs[m] );
  //     Op.Zdscal<N>( A[m], d_Ys[m] );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   // {
  //   //   CuC *d_sum, *d_dD_sum;
  //   //   CUDA_CHECK(cudaMalloc(&d_sum, N*CD));
  //   //   CUDA_CHECK(cudaMalloc(&d_dD_sum, N*CD));
  //   //   CUDA_CHECK(cudaMemcpy(d_sum, d_eta, N*CD, D2D));
  //   //   for(int m=1; m<size; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_sum, 1.0, d_Zs[m], d_sum);

  //   //   coo( d_dD_sum, d_sum );

  //   //   Dummy.dot<N>( &inner, d_eta, d_dD_sum );
  //   //   res += real(inner);

  //   //   CUDA_CHECK(cudaFree(d_sum));
  //   //   CUDA_CHECK(cudaFree(d_dD_sum));
  //   // }
  //   // {
  //   //   CuC *d_XZm, *d_dD_Ym;
  //   //   CUDA_CHECK(cudaMalloc(&d_XZm, N*CD));
  //   //   // CUDA_CHECK(cudaMalloc(&d_dD_Ym, N*CD));

  //   //   for(int m=1; m<size; m++) {
  //   //     {
  //   //       MatPoly Op;
  //   //       Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   //       Op.on_gpu<N>(d_XZm, d_Zs[m]);
  //   //     }
  //   //     coo( d_dD_Ym, d_Ys[m] );

  //   //     Dummy.dot<N>( &inner, d_XZm, d_dD_Ym );
  //   //     res -= real(inner);
  //   //   }

  //   //   CUDA_CHECK(cudaFree(d_XZm));
  //   //   CUDA_CHECK(cudaFree(d_dD_Ym));
  //   // }
  //   {
  //     CuC *d_XYm, *d_dD_Zm;
  //     CUDA_CHECK(cudaMalloc(&d_XYm, N*CD));
  //     CUDA_CHECK(cudaMalloc(&d_dD_Zm, N*CD));

  //     for(int m=1; m<size; m++) {
  //       {
  //         MatPoly Op;
  //         Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //         Op.on_gpu<N>(d_XYm, d_Ys[m]);
  //       }
  //       coo( d_dD_Zm, d_Zs[m] );

  //       Dummy.dot<N>( &inner, d_XYm, d_dD_Zm );
  //       res -= real(inner);
  //     }
  //     CUDA_CHECK(cudaFree(d_XYm));
  //     CUDA_CHECK(cudaFree(d_dD_Zm));
  //   }

  //   res *= 2.0*C/lambda_max;

  //   CUDA_CHECK(cudaFree(d_Xdag_eta));
  //   for(int m=1; m<size; m++) {
  //     CUDA_CHECK(cudaFree(d_Zs[m]));
  //     CUDA_CHECK(cudaFree(d_Ys[m]));
  //   }

  //   return res;
  // }


  // void mult_device2(CuC* d_res, const CuC* d_xi) const {
  //   const int mmax = size; // size;
  //   std::vector<CuC*> d_Zs(mmax);
  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));

  //   for(int m=1; m<mmax; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Zs[m], d_xi );
  //     Op.Zdscal<N>( A[m], d_Zs[m] );
  //   }

  //   // CUDA_CHECK(cudaMemcpy(d_res, d_xi, N*CD, D2D));
  //   CUDA_CHECK(cudaMemset(d_res, 0, N*CD));
  //   for(int m=1; m<mmax; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_Zs[m], d_res);

  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaFree(d_Zs[m]));
  // }



  // double grad2( const Link& link, const Gauge& U, CuC* d_xi ) const {
  //   const int mmax = size; // size
  //   std::vector<CuC*> d_Zs(mmax);
  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));

  //   for(int m=1; m<mmax; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Zs[m], d_xi );
  //   }

  //   COO coo;
  //   DW.d_coo_format(coo.en, U, link);
  //   coo.do_it();

  //   CuC *d_tmp1, *d_tmp2;
  //   CUDA_CHECK(cudaMalloc(&d_tmp1, N*CD));
  //   CUDA_CHECK(cudaMalloc(&d_tmp2, N*CD));

  //   double res = 0.0;
  //   MatPoly Dummy;

  //   for(int m=1; m<mmax; m++){
  //     coo( d_tmp1, d_Zs[m] ); // DH
  //     M_DWH( d_tmp2, d_tmp1 );

  //     CuC inner;
  //     Dummy.Zdscal<N>( A[m], d_Zs[m] );
  //     Dummy.dot<N>( &inner, d_Zs[m], d_tmp2 );
  //     res += -2.0 * real(inner) / (lambda_max * lambda_max);
  //   }

  //   CUDA_CHECK(cudaFree(d_tmp1));
  //   CUDA_CHECK(cudaFree(d_tmp2));
  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaFree(d_Zs[m]));

  //   return res;
  // }



// struct OverlapPseudoFermion {
//   using Complex = std::complex<double>;
//   using Link = std::array<int,2>; // <int,int>;
//   using Face = std::vector<int>;
//   using Gauge = U1onS2;
//   using Force=U1onS2;
//   using Rng = ParallelRng;
//   using Idx = long int;

//   static constexpr Complex I = Complex(0.0, 1.0);
//   const int NS=2;

//   const Dirac1fonS2& D;
//   const CGCUDA cg;
//   const Sparse& sparse;
//   Zolotarev sgn;

//   std::vector<Complex> phi;
//   std::vector<Complex> eta;

//   OverlapPseudoFermion()=delete;

//   OverlapPseudoFermion(const Dirac1fonS2& D_, const double k_=0.01, const int n_=21)
//     : D(D_)
//     , cg(D)
//     , sparse(cg.sparse)
//     , sgn(k_, n_)
//     , phi(D.lattice.n_sites*NS, 0.0)
//     , eta(D.lattice.n_sites*NS, 0.0)
//   {}

//   Complex operator[](const int i) const { return phi[i]; }
//   Complex& operator[](const int i) { return phi[i]; }

//   Complex operator()(const int ix, const int a) const { return phi[NS*ix+a]; }
//   Complex& operator()(const int ix, const int a) { return phi[NS*ix+a]; }

//   void H( Complex* res, const Complex* v, const U1onS2& U,
// 	  const double lambda_max = 1.0) const {
//     for(Idx i=0; i<sparse.N; i++) res[i] = v[i];

//     for(int m=1; m<sgn.size; m++) {
//       std::vector<Complex> tmp(sparse.N);
//       // cg.general_op_inv( tmp.data(), v, U, 1.0/(lambda_max*lambda_max), sgn.k*sgn.k/sgn.cp[m] );
//       cg.generalH_op_inv( tmp.data(), v, U, 1.0, sgn.k*sgn.k/sgn.cp[m], lambda_max );
//       for(Idx i=0; i<sparse.N; i++) res[i] += sgn.A[m]*tmp[i];
//     }

//     std::vector<Complex> tmp(sparse.N);
//     std::vector<Complex> tmp2(sparse.N);
//     for(Idx i=0; i<sparse.N; i++) tmp[i] = res[i];
//     multHW( tmp2, tmp, U, lambda_max );
//     const double aa = sgn.E * 2.0 / (1.0+sgn.lambda_inv) / (sgn.k*sgn.M);
//     for(Idx i=0; i<sparse.N; i++) res[i] = aa*tmp2[i];
//   }


//   void multD( std::vector<Complex>& Dxi, const std::vector<Complex>& xi, const Gauge& U ) const {
//     assert( Dxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csr[sparse.len];
//     D.coo_format(D_coo, U);
//     sparse.coo2csr( D_csr, D_coo );
//     sparse.mult<Complex>( Dxi.data(), xi.data(), D_csr );
//   }

//   void multHW( std::vector<Complex>& Dxi, const std::vector<Complex>& xi, const Gauge& U,
// 	       const double lambda_max ) const {
//     assert( Dxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csr[sparse.len];
//     D.H_coo_format(D_coo, U, lambda_max);
//     sparse.coo2csr( D_csr, D_coo );
//     sparse.mult<Complex>( Dxi.data(), xi.data(), D_csr );
//   }


//   void multDH( std::vector<Complex>& DHxi, const std::vector<Complex>& xi, const Gauge& U ) const {
//     assert( DHxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csrH[sparse.len];
//     D.coo_format(D_coo, U);
//     sparse.coo2csrH( D_csrH, D_coo );
//     sparse.multT<Complex>( DHxi.data(), xi.data(), D_csrH );    
//   }


//   void multDHD( std::vector<Complex>& DHDxi, const std::vector<Complex>& xi, const Gauge& U ) const {
//     assert( DHDxi.size()==D.lattice.n_sites*NS );
//     assert( xi.size()==D.lattice.n_sites*NS );

//     Complex D_coo[sparse.len], D_csr[sparse.len], D_csrH[sparse.len];
//     D.coo_format(D_coo, U);
//     sparse.coo2csr_csrH( D_csr, D_csrH, D_coo );

//     std::vector<Complex> tmp( D.lattice.n_sites*NS );
//     sparse.mult<Complex>( tmp.data(), xi.data(), D_csr );
//     sparse.multT<Complex>( DHDxi.data(), tmp.data(), D_csrH );    
//   }


// //   void gen( const Gauge& U, Rng& rng ) {
// //     const int N = D.lattice.n_sites*NS;
// //     std::vector<Complex> xi(N, 0.0);

// //     for(int ix=0; ix<D.lattice.n_sites; ix++) for(int a=0; a<NS; a++) xi[NS*ix+a] = ( rng.gaussian_site(ix) + I*rng.gaussian_site(ix) ) / std::sqrt(2.0);

// //     multDH( phi, xi, U );

// //     update_eta(U);
// //   }


// //   void update_eta( const Gauge& U ) { cg( eta.data(), phi.data(), U ); }


// //   auto begin(){ return phi.begin(); }
// //   auto end(){ return phi.end(); }
// //   auto begin() const { return phi.begin(); }
// //   auto end() const { return phi.end(); }


// //   Complex dot( const std::vector<Complex>& eta1, const std::vector<Complex>& xi) const {
// //     assert( eta1.size()==xi.size() );
// //     Complex res = 0.0;
// //     for(int i=0; i<eta1.size(); i++) res += std::conj(eta1[i]) * xi[i];
// //     return res;
// //   }

// //   Complex dot( const std::vector<Complex>& eta1 ) const {
// //     return dot(this->phi, eta1);
// //   }

// //   double S() const { return dot( eta ).real(); }


// //   double get_force( const Gauge& U, const Link& ell ) const {
// //     const int N = D.lattice.n_sites*NS;

// //     std::vector<Complex> dD;
// //     std::vector<int> is;
// //     std::vector<int> js;
// //     D.d_coo_format( dD, is, js, U, ell );

// //     std::vector<Complex> dD_eta(N);
// //     cg.sparse.multcoo( dD_eta, eta, dD, is, js );

// //     std::vector<Complex> DH_dD_eta(N);
// //     multDH( DH_dD_eta, dD_eta, U );

// //     return -2.0 * dot( eta, DH_dD_eta ).real();
// //   }


// //   Force dS( const Gauge& U ) const {
// //     Force pi( U.lattice ); // 0 initialized
// // #ifdef _OPENMP
// // #pragma omp parallel for num_threads(nparallel)
// // #endif
// //     for(int ell=0; ell<U.lattice.n_links; ell++) pi[ell] = get_force( U, U.lattice.links[ell] );
// //     return pi;
// //   }


// };







  // void mult_device3(CuC* d_res, const CuC* d_xi) const {
  //   const int mmax = size; // size;
  //   std::vector<CuC*> d_Zs(mmax);
  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaMalloc(&d_Zs[m], N*CD));

  //   for(int m=1; m<mmax; m++) {
  //     MatPoly Op;
  //     Op.push_back ( cplx(1.0/(lambda_max*lambda_max)), {&M_DW, &M_DWH} );

  //     const CuC a = cplx(-k*k/cp[m]);
  //     Op.push_back ( a, {} );
  //     Op.solve<N>( d_Zs[m], d_xi );
  //     Op.Zdscal<N>( A[m], d_Zs[m] );
  //   }

  //   CUDA_CHECK(cudaMemcpy(d_res, d_xi, N*CD, D2D));
  //   for(int m=1; m<mmax; m++) Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_Zs[m], d_res);

  //   CuC* d_tmp; CUDA_CHECK(cudaMalloc(&d_tmp, N*CD));
  //   MatPoly Op; Op.push_back ( cplx(1.0/(lambda_max)), {&M_DW} );
  //   CUDA_CHECK(cudaMemcpy(d_tmp, d_res, N*CD, D2D));
  //   Op.on_gpu<N>( d_res, d_tmp );

  //   Op.Zdscal<N>(C, d_res);
  //   Taxpy_gen<CuC,double,N><<<NBlocks, NThreadsPerBlock>>>(d_res, 1.0, d_xi, d_res); // 1+V

  //   for(int m=1; m<mmax; m++) CUDA_CHECK(cudaFree(d_Zs[m]));
  //   CUDA_CHECK(cudaFree(d_tmp));
  // }



  // auto f_update = std::bind(&Overlap::compute, &Dov, std::placeholders::_1);
  // LinOpWrapper M_mgrad_DHDov( f_mgrad_DHDov );
  // MatPoly Op_mgrad_DHDov;
  // Op_mgrad_DHDov.push_back ( cplx(1.0), {&M_DHDov} );

  CuC *d_eta, *d_xi;
  CUDA_CHECK(cudaMalloc(&d_eta, N*CD));
  CUDA_CHECK(cudaMalloc(&d_xi, N*CD));
  CUDA_CHECK(cudaMemcpy(d_xi, reinterpret_cast<const CuC*>(xi.data()), N*CD, H2D));

  {
    Dov.compute(U);
    std::cout << "# min max ratio: "
              << Dov.lambda_min << " "
              << Dov.lambda_max << " "
              << Dov.lambda_min/Dov.lambda_max << std::endl;
    std::cout << "# delta = " << Dov.Delta() << std::endl;

    Op_DHDov.solve<N>( d_eta, d_xi );
    Op_DHDov.dot<N>( &Sf, d_xi, d_eta );
    std::cout << "S = " << real(Sf) << " " << imag(Sf) << std::endl;
  }

  {
    Dov.compute(U);
    grad_d = f_mgrad_DHDov(ell, U, d_eta); // Dov.grad_device( ell, U, d_eta ); // DH
  }

  {
    Dov.compute(UP);
    Op_DHDov.solve<N>( d_eta, d_xi );
    Op_DHDov.dot<N>( &Sfp, d_xi, d_eta );
    std::cout << "Sp = " << real(Sfp) << " " << imag(Sfp) << std::endl;
  }
  {
    Dov.compute(UM);
    Op_DHDov.solve<N>( d_eta, d_xi );
    Op_DHDov.dot<N>( &Sfm, d_xi, d_eta );
    std::cout << "Sm = " << real(Sfm) << " " << imag(Sfm) << std::endl;
  }

  CuC check2 = (Sfp-Sfm)/(2.0*eps);
  std::cout << "grad = " << grad_d << std::endl;
  // std::cout << "grad = " << real(grad) << " " << imag(grad) << std::endl;
  std::cout << "check = " << real(check2) << " " << imag(check2) << std::endl;


  CUDA_CHECK(cudaFree(d_eta));
  CUDA_CHECK(cudaFree(d_xi));



  // auto f_Dov = std::bind(&Overlap::mult_device2, &Dov, std::placeholders::_1, std::placeholders::_2);
  // auto f_Dov = std::bind(&Overlap::mult_device3, &Dov, std::placeholders::_1, std::placeholders::_2);
  // auto f_Dov = std::bind(&Overlap::mult_device, &Dov, std::placeholders::_1, std::placeholders::_2);

  std::vector<Complex> Dxi(N), xi(N);
  for(int i=0; i<N; i++) xi[i] = rng.gaussian() + 1.0*Complex(0.0,1.0)*rng.gaussian();



  // ------------------

  const double M5 = -1.8;
  WilsonDirac DW(lattice, M5);

  Overlap Dov(DW, 11);

  const auto f_DHDov = std::bind(&Overlap::sq_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  LinOpWrapper M_DHDov( f_DHDov );
  MatPoly Op_DHDov; Op_DHDov.push_back ( cplx(1.0), {&M_DHDov} );
  auto f_DHov = std::bind(&Overlap::adj_device, &Dov, std::placeholders::_1, std::placeholders::_2);
  auto f_mgrad_DHDov = std::bind(&Overlap::grad_device, &Dov, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

  PseudoFermion pf( Op_DHDov, f_DHov, f_mgrad_DHDov );

  // ------------------

  Idx il=2;
  Link ell = lattice.links[il];

  const double eps = 1.0e-5;
  Gauge UP(U);
  UP[il] += eps;
  Gauge UM(U);
  UM[il] -= eps;

  Dov.compute(U);
  pf.gen( rng );
  double grad_d = pf.get_force( U, ell );

  Dov.compute(UP);
  pf.update_eta();
  double sfp = pf.S();

  Dov.compute(UM);
  pf.update_eta();
  double sfm = pf.S();

  double chck = (sfp-sfm)/(2.0*eps);
  std::cout << "grad = " << grad_d << std::endl;
  std::cout << "check = " << chck << std::endl;
