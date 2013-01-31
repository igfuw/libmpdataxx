/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_common.hpp"
#include "blitz.hpp"
#include "arakawa_c.hpp"

namespace solvers
{
  template<class bcx_t, class bcy_t, int n_eqs, typename real_t>
  class solver_2d : public solver_common<n_eqs, 2, real_t>
  {
    using parent = solver_common<n_eqs, 2, real_t>;
    using arr_2d_t = typename parent::arr_t;

    protected:
  
    bcx_t bcx;
    bcy_t bcy;

    int halo;
    rng_t i, j;

    void xchng(int e) // for current time level
    {
      bcx.fill_halos(this->psi[e][ this->n[e] ], j^halo);
      bcy.fill_halos(this->psi[e][ this->n[e] ], i^halo);
    }

    void xchng(int e, int lev) //for previous time level
    {
      bcx.fill_halos(this->psi[e][ this->n[e] - lev], j^halo);
      bcy.fill_halos(this->psi[e][ this->n[e] - lev], i^halo);
    }

    void xchng(arr_2d_t psi, rng_t range_i, rng_t range_j) // for a given array
    {
      bcx.fill_halos(psi, range_j);
      bcy.fill_halos(psi, range_i);
    }

    // ctor
    solver_2d(int nx, int ny, int halo) :
      halo(halo),
      i(0, nx-1), 
      j(0, ny-1),  
      bcx(rng_t(0, nx-1), halo), 
      bcy(rng_t(0, ny-1), halo)
    {
      for (int e = 0; e < n_eqs; ++e) // equations
        for (int n = 0; n < this->n_tlev; ++n) // time levels
          this->psi[e].push_back(new arr_2d_t(i^halo, j^halo));

      this->C.push_back(new arr_2d_t(i^h   , j^halo));
      this->C.push_back(new arr_2d_t(i^halo, j^h   ));
    }

    public:

    // accessor methods
    arr_2d_t state(int e = 0) 
    {
      return this->psi[e][ this->n[e] ](i,j).reindex({0,0});
    }
  };
}; // namespace solvers
