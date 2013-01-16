/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solvers_common.hpp"

namespace solvers
{
  template<class bcx_t, int n_eqs>
  class solver_1d : public solver_common<n_eqs>
  {
    //typedef arr_1d_t arr_t;

    protected:

    bcx_t bcx;
 
    // member fields
    arrvec_t<arr_1d_t> C;

    // psi contains model state including halo region
    arrvec_t<arr_1d_t> psi[n_eqs];

    int halo;
    rng_t i;

    void xchng() 
    {
      for (int e = 0; e < n_eqs; ++e) 
        bcx.fill_halos( psi[e][ this->n[e] ] );
    }

    // ctor
    solver_1d(int nx, int halo) :
      halo(halo),
      i(0, nx-1), 
      bcx(rng_t(0, nx-1), halo)
    {
      for (int e = 0; e < n_eqs; ++e) // equations
        for (int l = 0; l < 2; ++l) // time levels
          psi[e].push_back(new arr_1d_t(i^halo));

      C.push_back(new arr_1d_t(i^h));
    }

    public:

    // integration logic
    void solve(const int nt) 
    {
      for (int t = 0; t < nt; ++t) 
      {
        xchng();
        this->adv();
        this->cycle();
      }
    }

    // accessor method for psi (hides the halo region)
    arr_1d_t state(int e = 0) 
    {
      return psi[e][ this->n[e] ](i).reindex({0});
    }

    // accessor method for the Courant number field
    arr_1d_t courant() 
    { 
      return C[0]; 
    }
  };
}; // namespace solvers
