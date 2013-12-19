/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/solvers/adv/detail/solver_1d.hpp>
#include <libmpdata++/formulae/leapfrog_formulae.hpp>

// TODO: this doesn't work yet - just a draft! needs an initial non-leapfrog step

namespace libmpdataxx
{
  namespace solvers
  {
    template<
      typename real_t, 
      formulae::opts::opts_t opts,
      int halo = formulae::leapfrog::halo
    >
    class leapfrog_1d : public detail::solver<
      real_t, 
      1, 
      formulae::leapfrog::n_tlev, 
      opts,
      detail::max(halo, formulae::leapfrog::halo)
    > 
    {
      using parent_t = detail::solver<
        real_t, 
        1,  
        formulae::leapfrog::n_tlev, 
        opts,
        detail::max(halo, formulae::leapfrog::halo)
      >;

      void advop(int e)
      {
        formulae::leapfrog::op_1d( // TODO! G
          this->mem->psi[e], this->n[e], this->mem->GC[0], this->i
        );
      }

      public:

      struct params_t {};

      // ctor
      leapfrog_1d(
        typename parent_t::ctor_args_t args,
        const params_t &
      ) :
        parent_t(args) //mem, bcxl, bcxr, i)
      {}  
    };
  }; // namespace solvers
}; // namespace libmpdataxx

