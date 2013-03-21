/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <advoocat/solvers/detail/solver_1d.hpp>
#include <advoocat/formulae/leapfrog_formulae.hpp>

// TODO: this doesn't work yet - just a draft! needs an initial non-leapfrog step

namespace advoocat
{
  namespace solvers
  {
    template<
      typename real_t, 
      int n_eqs = 1, 
      int halo = formulae::leapfrog::halo
    >
    class leapfrog_1d : public detail::solver_1d<
      real_t, 
      1, 
      n_eqs, 
      formulae::leapfrog::n_tlev, 
      detail::max(halo, formulae::leapfrog::halo)
    > 
    {
      using parent_t = detail::solver_1d<
        real_t, 
        1,  
        n_eqs, 
        formulae::leapfrog::n_tlev, 
        detail::max(halo, formulae::leapfrog::halo)
      >;

      void advop(int e)
      {
        formulae::leapfrog::op_1d(
          this->mem->psi[e], this->n[e], this->mem->C[0], this->i
        );
      }

      public:

      struct params_t {};

      // ctor
      leapfrog_1d(
        typename parent_t::ctor_args_t args,
/*
        typename parent_t::mem_t *mem, 
        typename parent_t::bc_p &bcxl, 
        typename parent_t::bc_p &bcxr, 
        const rng_t &i, 
*/
        const params_t &
      ) :
        parent_t(args) //mem, bcxl, bcxr, i)
      {}  
    };
  }; // namespace solvers
}; // namespace advoocat

