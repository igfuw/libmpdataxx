/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_common.hpp"
#include "../blitz.hpp"
#include "../arakawa_c.hpp"

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      using namespace advoocat::arakawa_c;

      template<class bcx_t, class bcy_t, class mem_t, int n_tlev>
      class solver_2d : public solver_common<mem_t, n_tlev>
      {
	using parent_t = solver_common<mem_t, n_tlev>;
	using arr_2d_t = typename mem_t::arr_t;

	protected:
      
	bcx_t bcx;
	bcy_t bcy;

	int halo;
	rng_t i, j;

	void xchng(int e, int lev = 0) // for previous time levels
	{
          this->xchng(this->mem.psi[e][ this->mem.n[e] - lev], i^halo, j^halo);
	}

	void xchng(arr_2d_t psi, rng_t range_i, rng_t range_j) // for a given array
	{
	  bcx.fill_halos(psi, range_j);
	  bcy.fill_halos(psi, range_i);
	}

	// ctor
	solver_2d(mem_t &mem, const rng_t &i, const rng_t &j, const int halo) :
	  parent_t(mem),
	  halo(halo),
	  i(i), 
	  j(j),  
	  bcx(i, halo), 
	  bcy(j, halo)
	{}

	public:

	// empty by default
	static void alloc(mem_t &mem, const int nx, const int ny) 
        {
          const rng_t i(0, nx-1), j(0, ny-1);
          const int hlo = 1; //TODO

	  for (int e = 0; e < mem_t::n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem.psi[e].push_back(new arr_2d_t(i^hlo, j^hlo));

	  mem.C.push_back(new arr_2d_t(i^h   , j^hlo));
	  mem.C.push_back(new arr_2d_t(i^hlo, j^h   ));
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
