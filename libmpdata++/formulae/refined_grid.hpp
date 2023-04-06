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
      // averaged only from points within the domain, halos outside of the modeled domain are not used
      // however, halos between MPI domains are used;
      // it is assumed that these are halos are large enough!
      // currently, the halo between MPI domains is n_ref+n_ref/2
      template<class real_t, class arr_t, class ijk_t>
      void spatial_average_ref2reg(
        arr_t arr,
        const ijk_t &ijk,                    // index of the refined point that overlaps with the regular point
        const int &dist,                     // average calculate from the box (i-dist, i+dist)(j-dist,j+dist)(k-dist.k+dist)
        const std::array<int, 3> grid_size,  // number of cells in the domain in each direction (total, not for this MPI process), needed for averaging on edges so that it does not go out of the domain
        const bool volume_avg                // should averaging account weights proportional to cell volume (refined cells at box edges are only partially within the regular cell)
      )
      {
        // WARNING: very bad LOOPS
        for(int i = ijk[0].first(); i <= ijk[0].last(); i+= ijk[0].stride())
          for(int j = ijk[1].first(); j <= ijk[1].last(); j+= ijk[1].stride())
            for(int k = ijk[2].first(); k <= ijk[2].last(); k+= ijk[2].stride())
            {
              real_t norm = 0;
              for(int id=-dist; id<=dist; ++id)
                for(int jd=-dist; jd<=dist; ++jd)
                  for(int kd=-dist; kd<=dist; ++kd)
                  {
                    if(id==0 && jd==0 && kd==0) continue; // dont add the overlapping point
                    // dont add points outside of the domain
                    if(i+id < 0 || i+id > grid_size[0]-1 ||
                       j+jd < 0 || j+jd > grid_size[1]-1 ||
                       k+kd < 0 || k+kd > grid_size[2]-1) 
                       continue;

                    real_t w = 1; // weight
                    if(volume_avg)
                    {
                      // refined cells that are on the edge of the domain
                      if(i+id == 0 || i+id == grid_size[0]-1) w/=2.;
                      if(j+jd == 0 || j+jd == grid_size[1]-1) w/=2.;
                      if(k+kd == 0 || k+kd == grid_size[2]-1) w/=2.;

                      // refined cells that are half-between two resolved cells
                      if(id==-dist || id == dist) w/=2.;
                      if(jd==-dist || jd == dist) w/=2.;
                      if(kd==-dist || kd == dist) w/=2.;
                    }
                    arr(i,j,k) += w * arr(i+id,j+jd,k+kd);
                    norm += w;
                  }
              arr(i,j,k) /= norm+1; // +1 because the center cell was not added to norm
            }
      }
    };
  };
};
