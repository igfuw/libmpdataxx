/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_common.hpp"

namespace solvers
{
  template<class bcx_t, class bcy_t, class bcz_t, int n_eqs, typename real_t>
  class solver_3d : public solver_common<n_eqs, real_t>
  {
    protected:
  
    bcx_t bcx;
    bcy_t bcy;
    bcz_t bcz;

    // <TODO> this could in principle be placed in solver_common
    arrvec_t<arr_3d_t<real_t>> C, psi[n_eqs]; 
    int halo;
    // </TODO>

    rng_t i, j, k; // TODO: if stored as idx_t this also could be placed in solver_common

    void xchng(int e) 
    {
      bcx.fill_halos(psi[e][ this->n[e] ], j^halo, k^halo);
      bcy.fill_halos(psi[e][ this->n[e] ], k^halo, i^halo);
      bcz.fill_halos(psi[e][ this->n[e] ], i^halo, j^halo);
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
        for (int l = 0; l < 2; ++l) // time levels
          psi[e].push_back(new arr_3d_t<real_t>(i^halo, j^halo, k^halo));

      C.push_back(new arr_3d_t<real_t>(i^h, j^halo, k^halo));
      C.push_back(new arr_3d_t<real_t>(i^halo, j^h, k^halo));
      C.push_back(new arr_3d_t<real_t>(i^halo, j^halo, k^h));
    }

    public:

    // accessor methods
    arr_3d_t<real_t> state(int e = 0) 
    {
      return psi[e][ this->n[e] ](i,j,k).reindex({0,0,0});
    }

    arr_3d_t<real_t> courant(int d) // TODO: if arr_t would be a typedef in solve_common, this could be placed there as well
    { 
      return C[d]; 
    }
  };
}; // namespace solvers
