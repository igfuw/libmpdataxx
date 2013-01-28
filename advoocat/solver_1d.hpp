/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "solver_common.hpp"

namespace solvers
{
  template<class bcx_t, int n_eqs, typename real_t>
  class solver_1d : public solver_common<n_eqs, 1, real_t>
  {
    using parent = solver_common<n_eqs, 1, real_t>;
    using arr_1d_t = typename parent::arr_t;

    protected:

    bcx_t bcx;
 
    int halo;
    rng_t i;

    void xchng(int e) 
    {
      bcx.fill_halos( this->psi[e][ this->n[e] ] );
    }

    // ctor
    solver_1d(int nx, int halo) :
      halo(halo),
      i(0, nx-1), 
      bcx(rng_t(0, nx-1), halo)
    {
      for (int e = 0; e < n_eqs; ++e) // equations
        for (int n = 0; n < this->n_tlev; ++n) // time levels
          this->psi[e].push_back(new arr_1d_t(i^halo));

      this->C.push_back(new arr_1d_t(i^h));
    }

    public:

    // accessor method for psi (hides the halo region)
    arr_1d_t state(int e = 0) 
    {
      return this->psi[e][ this->n[e] ](i).reindex({0});
    }
  };
}; // namespace solvers
