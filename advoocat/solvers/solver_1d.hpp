/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "solver_common.hpp"

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      using namespace advoocat::arakawa_c;

      template<class bcx_t, class mem_t, int n_tlev>
      class solver_1d : public solver_common<mem_t, n_tlev>
      {
	using parent_t = solver_common<mem_t, n_tlev>;
	using arr_1d_t = typename mem_t::arr_t;

	protected:

	bcx_t bcx;
     
	int halo;
	rng_t i;

	void xchng(int e, int lev = 0) 
	{
	  bcx.fill_halos( this->mem.psi[e][ this->mem.n[e] - lev ] );
	}

	// ctor
	solver_1d(mem_t &mem, const rng_t &i, int halo) :
	  parent_t(mem), halo(halo), i(i), bcx(i, halo)
	{
	}

	public:

	static void alloc(mem_t &mem, const int nx)   
        {
          const rng_t i(0, nx-1);
          const int hlo = 1; // TODO!
          
	  for (int e = 0; e < mem_t::n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem.psi[e].push_back(new arr_1d_t(i^hlo));
    
	  mem.C.push_back(new arr_1d_t(i^h));
        } 
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
