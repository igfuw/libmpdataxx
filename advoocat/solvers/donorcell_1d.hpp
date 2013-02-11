/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "detail/solver_1d.hpp"
#include "../formulae/donorcell_formulae.hpp"

namespace advoocat
{
  namespace solvers
  {
    template<class bcx_t, class mem_t>
    class donorcell_1d : public detail::solver_1d<
      bcx_t, 
      mem_t, 
      formulae::donorcell::n_tlev, 
      formulae::donorcell::halo
    > 
    {
      using parent_t = detail::solver_1d<
        bcx_t, 
        mem_t, 
        formulae::donorcell::n_tlev, 
        formulae::donorcell::halo
      >;
   
      void advop(int e)
      {
        formulae::donorcell::op_1d(
          this->mem.psi[e], this->mem.n[e], this->mem.C[0], this->i
        );
      }

      public:

      struct params_t {};

      // ctor
      donorcell_1d(mem_t &mem, const rng_t &i, const params_t &) :
        parent_t(mem, i)
      { }  
    };
  }; // namespace solvers
}; // namespace advoocat

