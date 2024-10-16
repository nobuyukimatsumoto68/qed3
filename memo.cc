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
