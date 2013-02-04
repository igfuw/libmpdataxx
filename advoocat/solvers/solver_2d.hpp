/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_common.hpp"
#include "../blitz.hpp"
#include "../arakawa_c.hpp"

namespace advoocat
{
  namespace solvers
  {
    using namespace advoocat::arakawa_c;

// TODO: indent
  template<class bcx_t, class bcy_t, class mem_t>
  class solver_2d : public solver_common<mem_t>
  {
    using parent_t = solver_common<mem_t>;
    using arr_2d_t = typename mem_t::arr_t;

    protected:
  
    bcx_t bcx;
    bcy_t bcy;

    int halo;
    rng_t i, j;

    void xchng(int e) // for current time level
    {
      bcx.fill_halos(this->mem.psi[e][ this->mem.n[e] ], j^halo);
      bcy.fill_halos(this->mem.psi[e][ this->mem.n[e] ], i^halo);
    }

    void xchng(int e, int lev) //for previous time level
    {
      bcx.fill_halos(this->mem.psi[e][ this->mem.n[e] - lev], j^halo);
      bcy.fill_halos(this->mem.psi[e][ this->mem.n[e] - lev], i^halo);
    }

    void xchng(arr_2d_t psi, rng_t range_i, rng_t range_j) // for a given array
    {
      bcx.fill_halos(psi, range_j);
      bcy.fill_halos(psi, range_i);
    }

    // ctor
    solver_2d(mem_t &mem, const rng_t &i, const rng_t &j, const int halo) :
      parent_t(mem),
      halo(halo),
      i(i), 
      j(j),  
      bcx(i, halo), 
      bcy(j, halo)
    {
      for (int e = 0; e < mem_t::n_eqs; ++e) // equations
        for (int n = 0; n < this->n_tlev; ++n) // time levels
          mem.psi[e].push_back(new arr_2d_t(i^halo, j^halo));

      mem.C.push_back(new arr_2d_t(i^h   , j^halo));
      mem.C.push_back(new arr_2d_t(i^halo, j^h   ));
    }
    };
  }; // namespace solvers
}; // namespace advoocat
