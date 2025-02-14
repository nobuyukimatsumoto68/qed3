#pragma once

#include "s2.h"

struct S2Simp : public QfeLatticeS2 {
  using Link = std::array<Idx,2>; // <Idx,Idx>;
  using Face = std::vector<Idx>;

  std::vector<std::vector<Idx>> nns;

  std::map<const Link, const Idx> map2il;
  std::map<const Link, const int> map2sign;

  std::vector<Link> _links;

  std::vector<Idx> counter_accum;

  S2Simp(const int n_refine)
    : QfeLatticeS2(5, n_refine)
    , _links(n_links)
  {
    for(int ix=0; ix<this->n_sites; ix++){
      std::vector<Idx> tmp;
      for(int jj=0; jj<this->nn(ix); jj++) tmp.push_back( this->sites[ix].neighbors[jj] );
      nns.push_back(tmp);
    }
    {
      Idx counter=0;
      for(Idx ix=0; ix<this->n_sites; ix++){
        for(int jj=0; jj<this->nn(ix); jj++){
          counter++;
        }
        counter_accum.push_back(counter);
      }
    }

    {
      for(Idx ix=0; ix<this->n_sites; ix++){
        for(int jj=0; jj<this->nn(ix); jj++){
          const Idx iy = this->nns[ix][jj];
          if(iy>ix) continue;

          const int il = this->sites[ix].links[jj];
          _links[il] = Link{ix,iy}; // ix<iy
        }
      }

      for(Idx il=0; il<n_links; il++) {
	const Link link = _links[il];
	const Idx i = link[0];
	const Idx j = link[1];

	map2il.insert( { Link{i,j}, il } );
	map2il.insert( { Link{j,i}, il } );
      }

      for(Idx il=0; il<n_links; il++) {
	const Link link = _links[il];
	const Idx i = link[0];
	const Idx j = link[1];

	map2il.insert( { Link{i,j}, il } );
	map2il.insert( { Link{j,i}, il } );
      }

      for(Idx il=0; il<n_links; il++) {
	const auto link = _links[il];
	const Idx i = link[0];
	const Idx j = link[1];

	map2sign.insert( { Link{i,j}, +1 } );
	map2sign.insert( { Link{j,i}, -1 } );
      }
    }

  }

  inline int nn(const Idx ix) const {
    return this->sites[ix].nn;
  }
};
