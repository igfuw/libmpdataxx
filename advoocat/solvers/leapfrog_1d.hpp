/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "detail/solver_1d.hpp"
#include "../formulae/donorcell_formulae.hpp"

// TODO: this doesn't work yet - just a draft!

namespace advoocat
{
  namespace solvers
  {
    template<class bcx_t, class mem_t>
    class leapfrog_1d : public detail::solver_1d<bcx_t, mem_t, /* n_tlev = */ 3> 
    {
      void advop(int e)
      {
        formulae::leapfrog::op_1d(
          this->mem.psi[e], this->mem.n[e], this->mem.C[0], this->i
        );
      }

      public:

      struct params_t {};

      // ctor
      leapfrog_1d(mem_t &mem, const rng_t &i, const params_t &) :
        detail::solver_1d<bcx_t, mem_t, 3>(mem, i, /* halo = */ 1)
      { }  
    };
  }; // namespace solvers
}; // namespace advoocat

