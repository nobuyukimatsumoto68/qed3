  // isModdable
  std::cout << isModdable( 1.0e-15, 1.0 ) << std::endl;
  std::cout << isModdable( -1.0e-15, 1.0 ) << std::endl;
  std::cout << isModdable( 1.0+1.0e-15, 1.0 ) << std::endl;
  std::cout << isModdable( -1.0+1.0e-15, 1.0 ) << std::endl;
  std::cout << isModdable( 1.0-1.0e-15, 1.0 ) << std::endl;
  std::cout << isModdable( -1.0-1.0e-15, 1.0 ) << std::endl;
  std::cout << isModdable( 1.4-1.0e-15, 1.0 ) << std::endl;
  std::cout << isModdable( 1.4+1.0e-15, 1.0 ) << std::endl;

  // arccos
  std::cout << std::acos(1.0) << std::endl;
  std::cout << std::acos(-1.0) << std::endl;
  // std::cout << std::acos(1.0+1.0e-16) << std::endl;

  {
    // embedding, projection
    V3 x(0.0, 0.0, -2.0);
    V2 xi = projectionS2(x);
    V3 x_check = embedding3D(xi);
    std::cout << "x  = " << x.transpose() << std::endl;
    std::cout << "xi = " << xi.transpose() << std::endl;
    std::cout << "xc = " << x_check.transpose() << std::endl;
  }


  {
    // getSign
    I2 sign1, sign2;
    getSign(sign1, sign2, x1, x2);
    std::cout << sign1.transpose() << std::endl;
    std::cout << sign2.transpose() << std::endl;
  }

  // phi0
  {
    V2 xi1 = projectionS2(x1);
    V2 xi2 = projectionS2(x2);
    std::cout << "xi1 = " << xi1.transpose() << std::endl;
    std::cout << "xi2 = " << xi2.transpose() << std::endl;
    const double cosphi0 = getCosPhi0(xi1, xi2, -1, 1);
    double phi0 = getPhi0(xi1, xi2, -1, 1);
    std::cout << cosphi0 << std::endl;
    std::cout << phi0 << std::endl;
    bool is_sol_phi;
    getPhi0WithRatioCheck(phi0, is_sol_phi, xi1, xi2, 1, -1);
    std::cout << is_sol_phi << std::endl;
    std::cout << phi0 << std::endl;
  }



  // {
  //   Pt x1(0.130395, 0.990556, -0.042368);
  //   Pt x2(0.080886, 0.990556, 0.110921);

  //   std::cout << geodesicLength( x1, x2 ) << std::endl;
  // }

  // {
  //   Pt x1( 0.1, 1.0 );
  //   Pt x2( 0.5, 1.0 );

  //   Sol sol = SolveGeodesicsConstPhi( x1, x2 );

  //   std::cout << sol.theta(0.1) << std::endl;
  //   std::cout << sol.phi(0.1) << std::endl;
  // }

  // {
  //   Pt x1( 0.1, 0.1+M_PI );
  //   Pt x2( 0.5, 0.1 );

  //   Sol sol = SolveGeodesicsDeltaPhiEqPi( x1, x2 );

  //   std::cout << sol.theta(0.05, 0) << std::endl;
  //   std::cout << sol.phi(0.05, 0) << std::endl;

  //   std::cout << sol.theta(0.3, 1) << std::endl;
  //   std::cout << sol.phi(0.3, 1) << std::endl;
  // }

  // {
  //   Pt x1( 0.5*M_PI, 0.1 );
  //   Pt x2( 0.5, 0.1 );

  //   Sol sol = SolveGeodesicsEndPtPiHalf( x2, x1 );
  //   std::cout << sol.theta(0.1) << std::endl;
  //   std::cout << sol.phi(0.1) << std::endl;
  // }



  // {
  //   Pt x1( 0.5*M_PI, 0.1 );
  //   Pt x2( 0.5, 0.1 );

  //   Sol sol = SolveGeodesicsEndPtPiHalf( x2, x1 );
  //   std::cout << sol.theta(0.1) << std::endl;
  //   std::cout << sol.phi(0.1) << std::endl;
  // }

  // {
  //   Pt x1( 0.11594391015821674, 0.6448028608904725, 0.7555039909124803 );
  //   Pt x2( -0.00967100461303834, 0.5470756243846032, 0.8370273190726639 );

  //   Sol sol = SolveGeodesicsMonotonic( x1, x2 );
  //   std::cout << sol.theta(0.1) << std::endl;
  //   std::cout << sol.phi(0.1) << std::endl;
  // }

  // {
  //   Pt x1( -0.042367955582481, 0.13039515939997406, 0.990556438949753 );
  //   Pt x2( 0.11092074774878105, 0.08058864047764432, 0.990556438949753 );

  //   Sol sol = SolveGeodesicsAltering( x1, x2 );
  //   std::cout << sol.theta(0.05) << std::endl;
  //   std::cout << sol.phi(0.05) << std::endl;
  //   std::cout << sol.theta(0.15, 1) << std::endl;
  //   std::cout << sol.phi(0.15, 1) << std::endl;
  // }


  // Pt x1( -0.042367955582481, 0.13039515939997406, 0.990556438949753 );
  // Pt x2( 0.11092074774878105, 0.08058864047764432, 0.990556438949753 );
  Pt x1( 0.11594391015821674, 0.6448028608904725, 0.7555039909124803 );
  Pt x2( -0.00967100461303834, 0.5470756243846032, 0.8370273190726639 );

  Sol sol = SolveGeodesics( x1, x2 );

  std::cout << sol.theta( 0.5*sol.sE ) << std::endl;
  std::cout << sol.phi( 0.5*sol.sE ) << std::endl;

  if(sol.is_split){
    std::cout << sol.theta( 0.5*(sol.sE+sol.ell), 1) << std::endl;
    std::cout << sol.phi( 0.5*(sol.sE+sol.ell), 1) << std::endl;
  }



  // // monotonic
  // Pt x1( 0.11594391015821674, 0.6448028608904725, 0.7555039909124803 );
  // Pt x2( -0.00967100461303834, 0.5470756243846032, 0.8370273190726639 );
  // altering
  Pt x1( -0.042367955582481, 0.13039515939997406, 0.990556438949753 );
  Pt x2( 0.11092074774878105, 0.08058864047764432, 0.990556438949753 );

  Sol sol = SolveGeodesics( x1, x2 );

  std::cout << "theta = " << sol.theta( 0.5*sol.sE ) << std::endl;
  std::cout << "phi = " << sol.phi( 0.5*sol.sE ) << std::endl;

  std::cout << "Dtheta = "
	    << sol.Dtheta( 0.5*sol.sE ) << std::endl;
  std::cout << "Dtheta = "
	    << (sol.theta( 0.5*sol.sE+EPSNUMDER )-sol.theta( 0.5*sol.sE-EPSNUMDER ))/(2.0*EPSNUMDER) << std::endl;
  std::cout << "Dphi = "
	    << sol.Dphi( 0.5*sol.sE ) << std::endl;
  std::cout << "Dphi = "
	    << (sol.phi( 0.5*sol.sE+EPSNUMDER )-sol.phi( 0.5*sol.sE-EPSNUMDER ))/(2.0*EPSNUMDER) << std::endl;
  
  if(sol.is_split){
    double ss = 0.5*(sol.sE+sol.ell);
    std::cout << "theta = " << sol.theta( ss, 1) << std::endl;
    std::cout << "phi = " << sol.phi( ss, 1) << std::endl;

    std::cout << "Dtheta = "
	      << sol.Dtheta( ss, 1 ) << std::endl;
    std::cout << "Dtheta = "
	      << (sol.theta( ss+EPSNUMDER, 1)-sol.theta( ss-EPSNUMDER, 1))/(2.0*EPSNUMDER) << std::endl;
    std::cout << "Dphi = "
	      << sol.Dphi( ss ) << std::endl;
    std::cout << "Dphi = "
	      << (sol.phi( ss+EPSNUMDER, 1)-sol.phi( ss-EPSNUMDER, 1))/(2.0*EPSNUMDER) << std::endl;

  }
