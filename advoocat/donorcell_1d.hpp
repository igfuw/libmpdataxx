/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_1d.hpp"
#include "donorcell_common.hpp"

namespace solvers
{
  template<class bcx_t, int n_eqs = 1, typename real_t = float>
  class donorcell_1d : public solver_1d<bcx_t, n_eqs, real_t> 
  {
    void advop(int e)
    {
      donorcell::op_1d(
        this->psi[e], this->n[e], this->C[0], this->i
      );
    }

    public:

    // ctor
    donorcell_1d(int nx) :
      solver_1d<bcx_t, n_eqs, real_t>(nx, /* halo = */ 1)
    {}  
  };
}; // namespace solvers

