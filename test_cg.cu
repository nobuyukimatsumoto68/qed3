#include "rng.h"
#include "s2n.h"
#include "u1_s2_dual.h"
#include "dirac_s2_dual.h"
#include "sparse.h"
#include "cg_cuda.h"
#include "metropolis.h"


int main(int argc, char* argv[]){
  using MS=Eigen::Matrix2cd;
  using VD=Eigen::Vector2d;
  using VE=Eigen::Vector3d;
  using VC=Eigen::VectorXcd;
  using Complex=std::complex<double>;

  // geometry
  const int n_refine=2;

  Lattice lattice(n_refine);
  Dirac1fonS2 D(lattice);

  const double gR = 0.4;
  const double width = 5.0 * gR / std::sqrt( lattice.n_faces );

  const bool is_compact = false;
  U1onS2 U(lattice);
  U1Wilson SW(gR, is_compact);
  Metropolis<U1Wilson, U1onS2> met(SW, width);

  for(int i=0; i<100; i++) double r = met( U );

  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   // ----------------

  //   auto tmp1 = mat * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csr[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csr( v_csr, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];
  //   sparse.mult( res, v, v_csr );

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
  //   std::cout << "norm = " << norm << std::endl;
  // }


  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   // ----------------

  //   auto tmp1 = mat * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csr[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csr( v_csr, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];

  //   {
  //     const int N = sparse.N;
  //     const int len = sparse.len;

  //     CuC *d_val;
  //     cudacheck(cudaMalloc(&d_val, len*CD));
  //     cudacheck(cudaMemcpy(d_val, reinterpret_cast<CuC*>(v_csr), len*CD, H2D));
  //     //
  //     int *d_cols, *d_rows;
  //     cudacheck(cudaMalloc(&d_cols, len*sizeof(int)));
  //     cudacheck(cudaMalloc(&d_rows, (N+1)*sizeof(int)));
  //     cudacheck(cudaMemcpy(d_cols, sparse.cols_csr.data(), len*sizeof(int), H2D));
  //     cudacheck(cudaMemcpy(d_rows, sparse.rows_csr.data(), (N+1)*sizeof(int), H2D));
  //     //

  //     CuC *d_tmp, *d_v;
  //     cudacheck(cudaMalloc(&d_tmp, N*CD));
  //     cudacheck(cudaMemset(d_tmp, 0, N*CD));
  //     cudacheck(cudaMalloc(&d_v, N*CD));
  //     cudacheck(cudaMemcpy(d_v, reinterpret_cast<CuC*>(v), N*CD, H2D));

  //     mult<<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v,
  // 					  d_val, d_cols, d_rows, N);

  //     cudacheck(cudaMemcpy(res, reinterpret_cast<Complex*>(d_tmp), N*CD, D2H));

  //     cudacheck(cudaFree(d_val));
  //     cudacheck(cudaFree(d_cols));
  //     cudacheck(cudaFree(d_rows));
  //     cudacheck(cudaFree(d_tmp));
  //     cudacheck(cudaFree(d_v));
  //   }

  //   // sparse.mult( res, v, v_csr );

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++) {
  //     // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
  //     norm += std::abs(tmp1[i] - res[i]);
  //   }
  //   std::cout << "norm = " << norm << std::endl;

  // }


  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   // ----------------

  //   auto tmp1 = mat.adjoint() * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csrH[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csrH( v_csrH, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];

  //   {
  //     const int N = sparse.N;
  //     const int len = sparse.len;

  //     CuC *d_valH;
  //     cudacheck(cudaMalloc(&d_valH, len*CD));
  //     cudacheck(cudaMemcpy(d_valH, reinterpret_cast<CuC*>(v_csrH), len*CD, H2D));
  //     //
  //     int *d_colsT, *d_rowsT;
  //     cudacheck(cudaMalloc(&d_colsT, len*sizeof(int)));
  //     cudacheck(cudaMalloc(&d_rowsT, (N+1)*sizeof(int)));
  //     cudacheck(cudaMemcpy(d_colsT, sparse.cols_csrT.data(), len*sizeof(int), H2D));
  //     cudacheck(cudaMemcpy(d_rowsT, sparse.rows_csrT.data(), (N+1)*sizeof(int), H2D));
  //     //

  //     CuC *d_tmp, *d_v;
  //     cudacheck(cudaMalloc(&d_tmp, N*CD));
  //     cudacheck(cudaMemset(d_tmp, 0, N*CD));
  //     cudacheck(cudaMalloc(&d_v, N*CD));
  //     cudacheck(cudaMemcpy(d_v, reinterpret_cast<CuC*>(v), N*CD, H2D));

  //     mult<<<NBlocks, NThreadsPerBlock>>>(d_tmp, d_v,
  // 					  d_valH, d_colsT, d_rowsT, N);

  //     cudacheck(cudaMemcpy(res, reinterpret_cast<Complex*>(d_tmp), N*CD, D2H));

  //     cudacheck(cudaFree(d_valH));
  //     cudacheck(cudaFree(d_colsT));
  //     cudacheck(cudaFree(d_rowsT));
  //     cudacheck(cudaFree(d_tmp));
  //     cudacheck(cudaFree(d_v));
  //   }

  //   // sparse.mult( res, v, v_csr );

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++) {
  //     // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
  //     norm += std::abs(tmp1[i] - res[i]);
  //   }
  //   std::cout << "norm = " << norm << std::endl;

  // }

  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   auto tmp1 = (mat.adjoint()*mat).inverse() * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csr[ sparse.len ], v_csrH[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];

  //   solve(reinterpret_cast<CuC*>(res),
  // 	  reinterpret_cast<CuC*>(v),
  // 	  reinterpret_cast<CuC*>(v_csr), sparse.cols_csr, sparse.rows_csr,
  // 	  reinterpret_cast<CuC*>(v_csrH), sparse.cols_csrT, sparse.rows_csrT,
  // 	  sparse.N, sparse.len
  // 	  );

  //   // for(int i=0; i<sparse.N; i++){
  //   //   std::cout << res[i] << std::endl;
  //   // }

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++) {
  //     // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
  //     norm += std::abs(tmp1[i] - res[i]);
  //   }
  //   std::cout << "norm = " << norm << std::endl;

  // }



  {
    auto mat = D.matrix_form( U );
    const VC r = VC::Random( mat.cols() );
    auto tmp1 = (mat.adjoint()*mat).inverse() * r;

    // ----------------

    const CGCUDA solver( lattice, D );
    Complex v[solver.sparse.N];
    for(int i=0; i<solver.sparse.N; i++) v[i] = r[i];
    Complex res[solver.sparse.N];

    solver( res, v, U );

    double norm = 0.0;
    for(int i=0; i<solver.sparse.N; i++) {
      // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
      norm += std::abs(tmp1[i] - res[i]);
    }
    std::cout << "norm = " << norm << std::endl;

  }










  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   auto tmp1 = (mat.adjoint()*mat).inverse() * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csr[ sparse.len ], v_csrH[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];

  //   int *d_cols, *d_rows, *d_colsT, *d_rowsT;
  //   set_sparse( d_cols, d_rows, d_colsT, d_rowsT, sparse );

  //   solve(reinterpret_cast<CuC*>(res),
  // 	  reinterpret_cast<CuC*>(v),
  // 	  reinterpret_cast<CuC*>(v_csr), d_cols, d_rows,
  // 	  reinterpret_cast<CuC*>(v_csrH), d_colsT, d_rowsT,
  // 	  sparse.N, sparse.len
  // 	  );

  //   free_sparse( d_cols, d_rows, d_colsT, d_rowsT );

  //   // for(int i=0; i<sparse.N; i++){
  //   //   std::cout << res[i] << std::endl;
  //   // }

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++) {
  //     // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
  //     norm += std::abs(tmp1[i] - res[i]);
  //   }
  //   std::cout << "norm = " << norm << std::endl;

  // }



  // std::cout << "solve." << std::endl;
  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );
  //   auto tmp1 = (mat.adjoint()*mat).inverse() * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   CGCUDA cg(sparse);

  //   Complex v_coo[ sparse.len ], v_csr[ sparse.len ], v_csrH[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];

  //   cg.solve(reinterpret_cast<CuC*>(res),
  // 	     reinterpret_cast<CuC*>(v),
  // 	     reinterpret_cast<CuC*>(v_csr),
  // 	     reinterpret_cast<CuC*>(v_csrH)
  // 	     );

  //   // for(int i=0; i<sparse.N; i++){
  //   //   std::cout << res[i] << std::endl;
  //   // }

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++) {
  //     // std::cout << "i = " << i << ", " << std::abs(tmp1[i] - res[i]) << std::endl;
  //     norm += std::abs(tmp1[i] - res[i]);
  //   }
  //   std::cout << "norm = " << norm << std::endl;
  // }

  // std::cout << "solve end." << std::endl;

  


  
  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   // ----------------

  //   auto tmp1 = mat.transpose() * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csrT[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csrT( v_csrT, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];
  //   sparse.multT( res, v, v_csrT );

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
  //   std::cout << "norm = " << norm << std::endl;
  // }

  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   // ----------------

  //   auto tmp1 = mat.adjoint() * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csrH[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csrH( v_csrH, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];
  //   sparse.multT( res, v, v_csrH );

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
  //   std::cout << "norm = " << norm << std::endl;
  // }

  // {
  //   auto mat = D.matrix_form( U );
  //   const VC r = VC::Random( mat.cols() );

  //   // ----------------

  //   auto tmp1 = mat.adjoint() * mat * r;

  //   // ----------------

  //   Sparse sparse(lattice);
  //   Complex v_coo[ sparse.len ], v_csr[ sparse.len ], v_csrH[ sparse.len ];
  //   D.coo_format( v_coo, sparse.N, U );
  //   sparse.coo2csr_csrH( v_csr, v_csrH, v_coo );

  //   Complex res[sparse.N];
  //   Complex v[sparse.N];
  //   for(int i=0; i<sparse.N; i++) v[i] = r[i];
  //   sparse.mult( res, v, v_csr );
  //   for(int i=0; i<sparse.N; i++) v[i] = res[i];
  //   sparse.multT( res, v, v_csrH );

  //   double norm = 0.0;
  //   for(int i=0; i<sparse.N; i++)  norm += std::abs(tmp1[i] - res[i]);
  //   std::cout << "norm = " << norm << std::endl;
  // }

  return 0;
}
