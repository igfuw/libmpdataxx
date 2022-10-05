/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

//#include <libmpdata++/formulae/idxperm.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace refined
    {
      // calculate regular grid scalar as an average from the refined grid
      // modifies the refined grid scalar!
      // 3d version
      template<class real_t, class arr_t, class rng_t>
      void spatial_average_ref2reg(
        arr_t arr,
        const rng_t &i,  // index of the refined point that overlaps with the regular point
        const rng_t &j,  // ditto
        const rng_t &k,  // ditto
        const int &dist  // average calculate from the box (i-dist, i+dist)(j-dist,j+dist)(k-dist.k+dist)
      )
      {
//        using idxperm::pis;

        // debug output
         // std::cerr << "ranges in interpolation" << std::endl;
         // if(d==0)
         //   std::cerr << "range<" << d << ">: " << i << " " << j << " " << k << std::endl;
         // if(d==1)
         //   std::cerr << "range<" << d << ">: " << k << " " << i << " " << j << std::endl;
         // if(d==2)
         //   std::cerr << "range<" << d << ">: " << j << " " << k << " " << i << std::endl;

        for(int id=-dist; id<=dist; ++id)
          for(int jd=-dist; jd<=dist; ++jd)
            for(int kd=-dist; kd<=dist; ++kd)
            {
              if(id==0 && jd==0 && kd==0) continue; // dont add the overlapping point
              real_t w = 1; // weight of the point - points on edges of the averaging box have smaller weights
              if(id==-dist || id == dist) w/=2.;
              if(jd==-dist || jd == dist) w/=2.;
              if(kd==-dist || kd == dist) w/=2.;
              arr(i, j, k) += w * arr(i+id, j+jd, k+kd);
            }
        // divide by number of cells
        arr(i, j, k) /= pow(2*dist+1, 3);
      }
    };
  };
};
