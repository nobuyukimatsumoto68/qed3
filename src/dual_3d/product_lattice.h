#pragma once


template<typename BaseLattice>
struct ProductLattice {

  BaseLattice& base;
  const int Nt;

  ProductLattice( base& base_,
                  const int Nt_)
    : base(base_)
    , Nt(Nt_)
  {}

};
