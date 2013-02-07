/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_3d.hpp"
#include "../formulae/donorcell_formulae.hpp"

namespace advoocat
{
  namespace solvers
  {
    template<
      class bcx_t, 
      class bcy_t, 
      class bcz_t, 
      class mem_t
    > 
    class donorcell_3d : public detail::solver_3d<bcx_t, bcy_t, bcz_t, mem_t, /* n_tlev = */ 2> 
    {
      using parent_t = detail::solver_3d<bcx_t, bcy_t, bcz_t, mem_t, /* n_tlev = */ 2>;

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
	const rng_t &i, 
	const rng_t &j, 
	const rng_t &k, 
	const params_t &
      ) :
	parent_t(mem, i, j, k, /* halo = */ 1)
      {}  
    }; // class donorcell_3d
  }; // namespace solvers
}; // namespace advoocat
