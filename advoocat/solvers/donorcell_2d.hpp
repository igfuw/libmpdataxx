/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_2d.hpp"
#include "../formulae/donorcell_formulae.hpp"

namespace advoocat
{
  namespace solvers
  {
    template<class bcx_t, class bcy_t, class mem_t>
    class donorcell_2d : public detail::solver_2d<bcx_t, bcy_t, mem_t, /* n_tlev = */ 2> 
    {
      using parent_t = detail::solver_2d<bcx_t, bcy_t, mem_t, 2>;

      void advop(int e)
      {
        formulae::donorcell::op_2d(
          this->mem.psi[e], this->mem.n[e], this->mem.C, this->i, this->j
        );
      }

      public:

      struct params_t {};

      // ctor
      donorcell_2d(mem_t &mem, const rng_t &i, const rng_t &j, const params_t &) :
        parent_t(mem, i, j, /* halo = */ 1)
      {}  
    };
  }; // namespace solvers
}; // namespace advoocat

