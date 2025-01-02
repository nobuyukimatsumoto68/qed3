#pragma once

struct SparseMatrix{
  const Lattice& lattice;
  const int N;

  static constexpr int NS = 2;

  int len;

  std::vector<int> is;
  std::vector<int> js;

  std::vector<int> ell2em;
  std::vector<int> cols_csr;
  std::vector<int> rows_csr;

  std::vector<int> ell2emT;
  std::vector<int> cols_csrT;
  std::vector<int> rows_csrT;

  SparseHelper(const Lattice& lattice_)
    : lattice(lattice_)
    , N(NS*lattice.n_sites)
      };
