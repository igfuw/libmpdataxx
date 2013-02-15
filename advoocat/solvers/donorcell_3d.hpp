/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "detail/solver_3d.hpp"
#include "../formulae/donorcell_formulae.hpp"

namespace advoocat
{
  namespace solvers
  {
    template<
      typename real_t,
      int n_eqs = 1,
      int halo = formulae::donorcell::halo
    > 
    class donorcell_3d : public detail::solver_3d<
      real_t,
      3,
      n_eqs,
      formulae::donorcell::n_tlev, 
      detail::max(halo, formulae::donorcell::halo)
    >
    {
      using parent_t = detail::solver_3d<
        real_t,
        3,
        n_eqs,
        formulae::donorcell::n_tlev, 
        detail::max(halo, formulae::donorcell::halo)
      >;

      void advop(int e)
      {
	formulae::donorcell::op_3d(
	  this->mem->psi[e], this->n[e], this->mem->C, this->i, this->j, this->k
	);
      }

      public:

      struct params_t {};

      // ctor
      donorcell_3d(
	typename parent_t::mem_t *mem, 
        typename parent_t::bc_p &bcxl,
        typename parent_t::bc_p &bcxr,
        typename parent_t::bc_p &bcyl,
        typename parent_t::bc_p &bcyr,
        typename parent_t::bc_p &bczl,
        typename parent_t::bc_p &bczr,
	const rng_t &i, 
	const rng_t &j, 
	const rng_t &k, 
	const params_t &
      ) :
	parent_t(mem, bcxl, bcxr, bcyl, bcyr, bczl, bczr, i, j, k)
      {}  
    }; // class donorcell_3d
  }; // namespace solvers
}; // namespace advoocat
