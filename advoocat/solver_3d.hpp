/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_common.hpp"
#include "arakawa_c.hpp"

namespace solvers
{
  template<class bcx_t, class bcy_t, class bcz_t, int n_eqs, typename real_t>
  class solver_3d : public solver_common<n_eqs, 3, real_t>
  {
    using parent = solver_common<n_eqs, 3, real_t>;
    using arr_3d_t = typename parent::arr_t;

    protected:
  
    bcx_t bcx;
    bcy_t bcy;
    bcz_t bcz;

    // <TODO> this could in principle be placed in solver_common
    int halo;
    // </TODO>

    rng_t i, j, k; // TODO: if stored as idx_t this also could be placed in solver_common

    void xchng(int e) 
    {
      bcx.fill_halos(this->psi[e][ this->n[e] ], j^halo, k^halo);
      bcy.fill_halos(this->psi[e][ this->n[e] ], k^halo, i^halo);
      bcz.fill_halos(this->psi[e][ this->n[e] ], i^halo, j^halo);
    }

    // ctor
    solver_3d(int nx, int ny, int nz, int halo) :
      halo(halo),
      i(0, nx-1), 
      j(0, ny-1),  
      k(0, nz-1),  
      bcx(rng_t(0, nx-1), halo), 
      bcy(rng_t(0, ny-1), halo),
      bcz(rng_t(0, nz-1), halo)
    {
      for (int e = 0; e < n_eqs; ++e) // equations
        for (int n = 0; n < this->n_tlev; ++n) // time levels
          this->psi[e].push_back(new arr_3d_t(i^halo, j^halo, k^halo));

      this->C.push_back(new arr_3d_t(i^h, j^halo, k^halo));
      this->C.push_back(new arr_3d_t(i^halo, j^h, k^halo));
      this->C.push_back(new arr_3d_t(i^halo, j^halo, k^h));
    }

    public:

    // accessor methods
    arr_3d_t state(int e = 0) 
    {
      return this->psi[e][ this->n[e] ](i,j,k).reindex({0,0,0});
    }
  };
}; // namespace solvers
