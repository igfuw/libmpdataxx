/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_2d.hpp"
#include "donorcell_formulae.hpp"

namespace solvers
{
  template<class bcx_t, class bcy_t, int n_eqs = 1, typename real_t = float>
  class donorcell_2d : public solver_2d<bcx_t, bcy_t, n_eqs, real_t> 
  {
    void advop(int e)
    {
      donorcell::op_2d(
        this->psi[e], this->n[e], this->C, this->i, this->j
      );
    }

    public:

    // ctor
    donorcell_2d(int nx, int ny) :
      solver_2d<bcx_t, bcy_t, n_eqs, real_t>(nx, ny, /* halo = */ 1)
    {}  
  };
}; // namespace solvers

