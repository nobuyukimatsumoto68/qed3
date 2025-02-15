



  // void set_ell_and_link_volumes() { // geodesic
  //   for(int il=0; il<lattice.n_links; il++) {
  //     const auto link = lattice.links[il];
  //     const int iA = link.faces[0];
  //     const int iB = link.faces[1];

  //     double ellA=0.0, ellB=0.0;
  //     double areaA=0.0, areaB=0.0;
  //     {
  // 	const QfeFace& face = lattice.faces[iA];
  // 	Vec3 r0, r1, r2; // r0,1: link
  // 	if(face.sites[0]==link.sites[0] && face.sites[1]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[0]];
  // 	  r1 = lattice.r[face.sites[1]];
  // 	  r2 = lattice.r[face.sites[2]];
  // 	}
  // 	else if(face.sites[1]==link.sites[0] && face.sites[2]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[1]];
  // 	  r1 = lattice.r[face.sites[2]];
  // 	  r2 = lattice.r[face.sites[0]];
  // 	}
  // 	else if(face.sites[2]==link.sites[0] && face.sites[0]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[2]];
  // 	  r1 = lattice.r[face.sites[0]];
  // 	  r2 = lattice.r[face.sites[1]];
  // 	} // reverse
  // 	else if(face.sites[1]==link.sites[0] && face.sites[0]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[1]];
  // 	  r1 = lattice.r[face.sites[0]];
  // 	  r2 = lattice.r[face.sites[2]];
  // 	}
  // 	else if(face.sites[2]==link.sites[0] && face.sites[1]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[2]];
  // 	  r1 = lattice.r[face.sites[1]];
  // 	  r2 = lattice.r[face.sites[0]];
  // 	}
  // 	else if(face.sites[0]==link.sites[0] && face.sites[2]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[0]];
  // 	  r1 = lattice.r[face.sites[2]];
  // 	  r2 = lattice.r[face.sites[1]];
  // 	}
  // 	else assert(false);

  // 	//
  // 	const Vec3 p = circumcenter(r0, r1, r2).transpose();
  // 	assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
  // 	assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

  // 	double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
  // 	double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
  // 	double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

  // 	double s_ = 0.5*(a_+b_+c_);
  // 	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
  // 	double area_ = 4.0*std::atan( std::sqrt( tmp ) );

  // 	ellA = c_;
  // 	areaA = area_;
  //     }
  //     {
  // 	const QfeFace& face = lattice.faces[iB];
  // 	Vec3 r0, r1, r2; // r0,1: link
  // 	if(face.sites[0]==link.sites[0] && face.sites[1]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[0]];
  // 	  r1 = lattice.r[face.sites[1]];
  // 	  r2 = lattice.r[face.sites[2]];
  // 	}
  // 	else if(face.sites[1]==link.sites[0] && face.sites[2]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[1]];
  // 	  r1 = lattice.r[face.sites[2]];
  // 	  r2 = lattice.r[face.sites[0]];
  // 	}
  // 	else if(face.sites[2]==link.sites[0] && face.sites[0]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[2]];
  // 	  r1 = lattice.r[face.sites[0]];
  // 	  r2 = lattice.r[face.sites[1]];
  // 	} // reverse
  // 	else if(face.sites[1]==link.sites[0] && face.sites[0]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[1]];
  // 	  r1 = lattice.r[face.sites[0]];
  // 	  r2 = lattice.r[face.sites[2]];
  // 	}
  // 	else if(face.sites[2]==link.sites[0] && face.sites[1]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[2]];
  // 	  r1 = lattice.r[face.sites[1]];
  // 	  r2 = lattice.r[face.sites[0]];
  // 	}
  // 	else if(face.sites[0]==link.sites[0] && face.sites[2]==link.sites[1]){
  // 	  r0 = lattice.r[face.sites[0]];
  // 	  r1 = lattice.r[face.sites[2]];
  // 	  r2 = lattice.r[face.sites[1]];
  // 	}
  // 	else assert(false);

  // 	const Vec3 p = circumcenter(r0, r1, r2).transpose();
  // 	assert( std::abs( (p-r0).norm() - (p-r1).norm() )<1.0e-14 );
  // 	assert( std::abs( (p-r0).norm() - (p-r2).norm() )<1.0e-14 );

  // 	double a_ = std::acos( r0.dot(p) /(r0.norm()* p.norm()) );
  // 	double b_ = std::acos( r1.dot(p) /(r1.norm()* p.norm()) );
  // 	double c_ = std::acos( r0.dot(r1)/(r0.norm()*r1.norm()) ); // ell

  // 	double s_ = 0.5*(a_+b_+c_);
  // 	double tmp = std::tan(0.5*s_) * std::tan(0.5*(s_-a_)) * std::tan(0.5*(s_-b_)) * std::tan(0.5*(s_-c_));
  // 	double area_ = 4.0*std::atan( std::sqrt( tmp ) );

  // 	ellB = c_;
  // 	areaB = area_;
  //     }

  //     assert( std::abs(ellA-ellB)<1.0e-14 );
  //     ell[il] = ellA;
  //     link_volume[il] = areaA + areaB;
  //   }

  //   int counter = 0;
  //   a = 0.0;
  //   for(int il=0; il<lattice.n_links; il++) {
  //     a += ell[il];
  //     counter++;
  //     std::cout << "ell[il] =" << ell[il] << std::endl;
  //   }
  //   a /= counter;
  //   a *= 1.0;
  //   std::cout << "a = " << a << std::endl;
  // }
























  // void set_site_vol(){
  //   for(int i=0; i<lattice.n_sites; i++){
  //     site_vol[i] = 0.0;
  //     const auto x = lattice.sites[i];
  //     //for(const int il : x.links){
  //     for(int jj=0; jj<x.nn; jj++){
  // 	const int il = x.links[jj];
  // 	site_vol[i] += 0.25*ell[il]*ellstar[il];
  //     }
  //   }
  // }
  
