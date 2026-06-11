#pragma once

// sparse_dirac_claude.h
// -----------------------------------------------------------------------------
// Copy of sparse_dirac.h with ONE change in DWDevice::initialize(): the COO->CSR row inversion is
// done by an O(len) bucketing (counting sort) instead of the original O(N*len) double scan, which
// dominated startup ("SCR indexing" setup) at larger L (N*nnz ~ 5e10 at L=4).
//
// The bucketing is byte-identical to the original (ell stays ascending within each row).  The
// original scan is KEPT below under -DCSR_VERIFY: when defined, it is rebuilt into check arrays and
// asserted equal to the bucketing result (one-time verification / rollback reference).  A build-time
// print reports the index-build wall time (N, nnz) so the speedup is visible in the log.
//
// Same struct name (DWDevice) -- consumers swap `#include "sparse_dirac.h"` for this and recompile.
// After verification, switch the include everywhere.
// -----------------------------------------------------------------------------


template<class WilsonDirac>
struct DWDevice{
  using T = CuC;

  const WilsonDirac& D;
  const Idx N;

  static constexpr int NS = 2;

  Idx len;

  Idx *ell2em, *ell2emT, *cols_csr, *rows_csr, *cols_csrT, *rows_csrT;
  Complex *v_coo, *v_csr, *v_csrH;

  Idx *d_cols, *d_rows, *d_colsT, *d_rowsT;
  CuC *d_val, *d_valH;

  std::vector<Idx> is;
  std::vector<Idx> js;


  DWDevice(const WilsonDirac& D_)
    : D(D_)
    , N(Comp::N)
  {
    initialize();
  } // end of constructor


  ~DWDevice()
  {
    CUDA_CHECK(cudaFree(d_cols));
    CUDA_CHECK(cudaFree(d_rows));
    CUDA_CHECK(cudaFree(d_colsT));
    CUDA_CHECK(cudaFree(d_rowsT));

    CUDA_CHECK(cudaFree(d_val));
    CUDA_CHECK(cudaFree(d_valH));

    CUDA_CHECK(cudaFreeHost(ell2em));
    CUDA_CHECK(cudaFreeHost(ell2emT));
    CUDA_CHECK(cudaFreeHost(cols_csr));
    CUDA_CHECK(cudaFreeHost(rows_csr));
    CUDA_CHECK(cudaFreeHost(cols_csrT));
    CUDA_CHECK(cudaFreeHost(rows_csrT));
    CUDA_CHECK(cudaFreeHost(v_coo));
    CUDA_CHECK(cudaFreeHost(v_csr));
    CUDA_CHECK(cudaFreeHost(v_csrH));
  }

  // @@@
  void initialize(){
    // ========= COO ========= //

    D.coo_structure(is, js);
    len = js.size();
    assert( is.size()==len );

    // ========= CSR ========= //

    Idx em=0;
    Idx emT=0;

    CUDA_CHECK( cudaMallocHost( &ell2em, len*sizeof(Idx) ) );
    CUDA_CHECK( cudaMallocHost( &cols_csr, len*sizeof(Idx) ) );
    memset(ell2em, 0, len*sizeof(Idx));
    memset(cols_csr, 0, len*sizeof(Idx));

    CUDA_CHECK( cudaMallocHost( &rows_csr, (N+1)*sizeof(Idx) ) );
    memset(rows_csr, 0, (N+1)*sizeof(Idx));

    CUDA_CHECK( cudaMallocHost( &ell2emT, len*sizeof(Idx) ) );
    CUDA_CHECK( cudaMallocHost( &cols_csrT, len*sizeof(Idx) ) );
    memset(ell2emT, 0, len*sizeof(Idx));
    memset(cols_csrT, 0, len*sizeof(Idx));

    CUDA_CHECK( cudaMallocHost( &rows_csrT, (N+1)*sizeof(Idx) ) );
    memset(rows_csrT, 0, (N+1)*sizeof(Idx));

    Idx counter=0;
    rows_csr[counter] = em;
    rows_csrT[counter] = em;
    counter++;


    Timer csr_timer;   // _claude: benchmark the CSR index build (was the O(N*len) bottleneck)

    std::vector<std::vector<Idx>> ell_record(N);
    std::vector<std::vector<Idx>> ell_recordT(N);

    // _claude: O(len) bucketing of the COO entries by row (a counting sort).  For each entry ell,
    // append it to its source-row list (ell_record[is[ell]]) and column-row list (ell_recordT[js[ell]]).
    // Iterating ell in ascending order keeps each per-row list ascending in ell, so the lists are
    // byte-identical to the original O(N*len) double scan -- but ~N x faster (no full sweep per row).
    for(Idx ell=0; ell<len; ell++){
      ell_record[ is[ell] ].push_back(ell);
      ell_recordT[ js[ell] ].push_back(ell);
    }

    // ORIGINAL O(N*len) double scan (sparse_dirac.h), kept for rollback + verification.  With
    // -DCSR_VERIFY it is rebuilt into separate arrays and asserted identical to the bucketing above.
#ifdef CSR_VERIFY
    {
      std::vector<std::vector<Idx>> chk(N), chkT(N);
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL_SORT)
#endif
      for(Idx i=0; i<N; i++){
        for(Idx ell=0; ell<len; ell++){
          if( is[ell]==i ) chk[i].push_back(ell);
          if( js[ell]==i ) chkT[i].push_back(ell);
        }
      }
      for(Idx i=0; i<N; i++){ assert(chk[i]==ell_record[i]); assert(chkT[i]==ell_recordT[i]); }
      std::cout << "# CSR_VERIFY: O(len) bucketing == O(N*len) scan  (N=" << N << ", nnz=" << len << ")\n";
    }
#endif

    for(Idx i=0; i<N; i++){
      for(const Idx ell : ell_record[i]){
        ell2em[ell] = em;
        cols_csr[em] = js[ell];
        ++em;
      }
      for(const Idx ell : ell_recordT[i]){
        ell2emT[ell] = emT;
        cols_csrT[emT] = is[ell];
        ++emT;
      }
      rows_csr[counter] = em;
      rows_csrT[counter] = emT;
      counter++;
    }

    assert(counter==N+1);

    std::cout << "# CSR index build (_claude bucketing): " << csr_timer.currentSeconds()
              << " s  (N=" << N << ", nnz=" << len << ")\n";

    CUDA_CHECK(cudaMalloc(&d_cols, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rows, (N+1)*sizeof(Idx)));
    CUDA_CHECK(cudaMemcpy(d_cols, cols_csr, len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rows, rows_csr, (N+1)*sizeof(Idx), H2D));
    //
    CUDA_CHECK(cudaMalloc(&d_colsT, len*sizeof(Idx)));
    CUDA_CHECK(cudaMalloc(&d_rowsT, (N+1)*sizeof(Idx)));
    CUDA_CHECK(cudaMemcpy(d_colsT, cols_csrT, len*sizeof(Idx), H2D));
    CUDA_CHECK(cudaMemcpy(d_rowsT, rows_csrT, (N+1)*sizeof(Idx), H2D));
    //
    CUDA_CHECK(cudaMalloc(&d_val, len*CD));
    CUDA_CHECK(cudaMalloc(&d_valH, len*CD));

    CUDA_CHECK( cudaMallocHost( &v_coo, len*CD ) );
    CUDA_CHECK( cudaMallocHost( &v_csr, len*CD ) );
    CUDA_CHECK( cudaMallocHost( &v_csrH, len*CD ) );
    memset(v_coo, 0, len*CD);
    memset(v_csr, 0, len*CD);
    memset(v_csrH, 0, len*CD);
  }

  void coo2csr_csrH() const {
#ifdef _OPENMP
#pragma omp parallel for num_threads(Comp::NPARALLEL)
#endif
    for(Idx ell=0; ell<len; ell++) {
      v_csr[ ell2em[ell] ] = v_coo[ell];
      v_csrH[ ell2emT[ell] ] = std::conj( v_coo[ell] );
    }
  }

  template<typename Gauge>
  void update( const Gauge& U ){
    D.coo_format( v_coo, U );
    coo2csr_csrH();

    CUDA_CHECK(cudaMemcpy(d_val, v_csr, len*CD, H2D));
    CUDA_CHECK(cudaMemcpy(d_valH, v_csrH, len*CD, H2D));
  }


  void associateCSR( CSR<Comp::N>& M, const bool is_transpose=false ){
    if(!is_transpose){
      M.cols = this->d_cols;
      M.rows = this->d_rows;
      M.val = this->d_val;
    }
    else{
      M.cols = this->d_colsT;
      M.rows = this->d_rowsT;
      M.val = this->d_valH;
    }
  }
};
