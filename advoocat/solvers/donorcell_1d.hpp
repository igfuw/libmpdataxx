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
    template<class mem_t, int halo = formulae::donorcell::halo>
    class donorcell_1d : public detail::solver_1d<
      mem_t, 
      formulae::donorcell::n_tlev, 
      detail::max(halo, formulae::donorcell::halo)
    > 
    {
      using parent_t = detail::solver_1d<
        mem_t, 
        formulae::donorcell::n_tlev, 
        detail::max(halo, formulae::donorcell::halo)
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
      donorcell_1d(
        mem_t &mem, 
        typename parent_t::bc_p &bcxl, 
        typename parent_t::bc_p &bcxr, 
        const rng_t &i, 
        const params_t &
      ) :
        parent_t(mem, bcxl, bcxr, i)
      { }  
    };
  }; // namespace solvers
}; // namespace advoocat

