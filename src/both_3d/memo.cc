        // for(int il=0; il<base.n_links; il++){
        //   if( std::abs( base.ell[il]-base.ell[il0] )>1.0e-10 ) continue;
        //   // if( std::abs( base.ell[il]-base.ell[il0] )<1.0e-10 ) continue;
        //   counter2++;
        //   // tmp2[s] += std::cos( U.plaquette_angle(s, U.lattice.links[il]) );
        //   // std::cout << "debug. angle = " << U.plaquette_angle(s, U.lattice.links[il]) << std::endl;
        //   // std::cout << "debug. factor = " << std::pow(base.mean_ell/base.ell[il], 2) << std::endl;
        //   tmp2[s] += std::pow(base.mean_ell/base.ell[il], 2) * ( std::cos( U.plaquette_angle(s, U.lattice.links[il]) ) - 1.0);
        // }
        // tmp1[s] /= counter1;
        // tmp2[s] /= counter2;
        // std::cout << "debug. tmp2[s] = " << tmp2[s] << std::endl;
        // std::cout << "debug. counter2 = " << counter2 << std::endl;


            // tmp1[s] += std::cos( U.plaquette_angle(s, U.lattice.faces[i_face]) );
            // tmp1[s] += std::pow(base.mean_vol/base.vols[i_face], 2) * ( std::cos( U.plaquette_angle(s, U.lattice.faces[i_face]) ) - 1.0);


      // for(int s=1; s<Comp::Nt; s++){
      //   tmp1[0] += tmp1[s];
      //   tmp2[0] += tmp2[s];
      // }


      // for(int t=0; t<Comp::Nt; t++){
      //   tmp1[t]/=Comp::Nt;
      //   plaq_s0[t].push_back( tmp1[0]/Comp::Nt );
      // }

      // plaq_t0.push_back( tmp2[0]/Comp::Nt );
      // std::cout << "debug. " << tmp2[0]/Comp::Nt << std::endl;


  double mean_s0=0.0, var_s0=0.0;
  double mean_t0=0.0, var_t0=0.0;
  for(int k=0; k<plaq_s0.size(); k++){
    mean_s0 += plaq_s0[k];
    mean_t0 += plaq_t0[k];
  }
  mean_s0 /= plaq_s0.size();
  mean_t0 /= plaq_s0.size();

  for(int k=0; k<plaq_s0.size(); k++){
    var_s0 += (plaq_s0[k]-mean_s0)*(plaq_s0[k]-mean_s0);
    var_t0 += (plaq_t0[k]-mean_t0)*(plaq_t0[k]-mean_t0);
  }
  var_s0 /= plaq_s0.size()*plaq_s0.size();
  var_t0 /= plaq_t0.size()*plaq_t0.size();

  // std::cout << "s0: " << beta << ", " << mean_s0 << ", " << std::sqrt( var_s0 ) << ", " << 0.5 * base.mean_vol/base.vols[iface0] * SW.beta_s << ", " << mean_s0 - 0.5 * base.mean_vol/base.vols[iface0] * SW.beta_s << std::endl;
  // std::cout << "t0: " << beta << ", " << mean_t0 << ", " << std::sqrt( var_t0 ) << ", " << 0.5 *base.mean_ell/base.ell[il0]* SW.beta_t << ", " << mean_t0 - 0.5 *base.mean_ell/base.ell[il0]* SW.beta_t << std::endl;
  std::cout << "s0: " << iface0 << " "
            << beta << ", " << mean_s0 << ", " << std::sqrt( var_s0 ) << std::endl;
  std::cout << "t0: " << beta << ", " << mean_t0 << ", " << std::sqrt( var_t0 ) << std::endl;
  // std::cout << "factor = " << U.lattice.mean_vol/U.lattice.vols[iface0] << std::endl;
