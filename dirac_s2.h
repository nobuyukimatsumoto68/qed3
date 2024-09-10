#pragma once

#include "s2.h"

struct Dirac1fonS2 {
  QfeLatticeS2& lattice;

  Dirac1fonS2()=delete;

  Dirac1fonS2(QfeLatticeS2& lattice_)
    : lattice(lattice_)
  {}

  Dirac1fonS2( const Dirac1fonS2& other )
    : lattice(other.lattice)
  {}

  Dirac1fonS2 & operator=(const Dirac1fonS2&) = delete;


};
