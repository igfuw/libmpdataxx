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
      class mem_t,
      int halo = formulae::donorcell::halo
    > 
    class donorcell_3d : public detail::solver_3d<
      mem_t, 
      formulae::donorcell::n_tlev, 
      detail::max(halo, formulae::donorcell::halo)
    >
    {
      using parent_t = detail::solver_3d<
        mem_t, 
        formulae::donorcell::n_tlev, 
        detail::max(halo, formulae::donorcell::halo)
      >;

      void advop(int e)
      {
	formulae::donorcell::op_3d(
	  this->mem.psi[e], this->mem.n[e], this->mem.C, this->i, this->j, this->k
	);
      }

      public:

      struct params_t {};

      // ctor
      donorcell_3d(
	mem_t &mem, 
        typename parent_t::bc_p &bcx,
        typename parent_t::bc_p &bcy,
        typename parent_t::bc_p &bcz,
	const rng_t &i, 
	const rng_t &j, 
	const rng_t &k, 
	const params_t &
      ) :
	parent_t(mem, bcx, bcy, bcz, i, j, k)
      {}  
    }; // class donorcell_3d
  }; // namespace solvers
}; // namespace advoocat
