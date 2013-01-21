/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_3d.hpp"
#include "donorcell_formulae.hpp"

namespace solvers
{
  template<
    class bcx_t, 
    class bcy_t, 
    class bcz_t, 
    int n_eqs = 1, 
    typename real_t = float
  > 
  class donorcell_3d : public solver_3d<bcx_t, bcy_t, bcz_t, n_eqs, real_t> 
  {
    void advop(int e)
    {
      donorcell::op_3d(
        this->psi[e], this->n[e], this->C, this->i, this->j, this->k
      );
    }

    public:

    // ctor
    donorcell_3d(int nx, int ny, int nz) :
      solver_3d<bcx_t, bcy_t, bcz_t, n_eqs, real_t>(nx, ny, nz, /* halo = */ 1)
    {}  
  };
}; // namespace solvers

