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

  double alat;

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
        // Idx counter_tmp=0;
        counter_accum.push_back(counter);
        for(int jj=0; jj<this->nn(ix); jj++) counter+=8;
        // counter += counter_tmp*8;
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


    {
      double mean_vol=0.;
      Idx counter=0;
      std::cout << "# reading vols" << std::endl;
      std::ifstream file(dir+"dualtriangleareas_n"+std::to_string(n_refine)+"_singlepatch.dat");
      assert(file.is_open());
      std::string str;
      while (std::getline(file, str)){
	std::istringstream iss(str);
	double v1;
	iss >> v1;
	// vol.push_back(v1);
        mean_vol += v1;
	counter++;
      }
      mean_vol /= counter;
      alat = std::sqrt( mean_vol*4.0/std::sqrt(3.0) );
    }
    // alat = std::sqrt( 8.0*M_PI/std::sqrt(3.0)/n_sites );
    // alat = std::sqrt( 8.0*M_PI/std::sqrt(3.0)/this->n_faces );
  }

  inline int nn(const Idx ix) const {
    return this->sites[ix].nn;
  }
};
