#pragma once

template<typename T>
struct SparseMatrix{
  Idx len;

  std::vector<T> val;
  std::vector<Idx> cols;
  std::vector<Idx> rows;

  SparseMatrix(const Idx len_)
    : len(len_)
  {}

};


