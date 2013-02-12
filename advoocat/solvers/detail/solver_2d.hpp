/** @file 
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_common.hpp"
#include "../../arakawa_c.hpp"
#include "../../bcond/bcond.hpp"

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      using namespace advoocat::arakawa_c;

      template<class mem_t, int n_tlev, int halo>
      class solver_2d : public solver_common<mem_t, n_tlev, halo>
      {
	using parent_t = solver_common<mem_t, n_tlev, halo>;
	using arr_2d_t = typename mem_t::arr_t;

	protected:
      
        typedef std::unique_ptr<bcond::bcond_t<typename mem_t::real_t>> bc_p; // TODO: move to parent
	bc_p bcxl, bcxr, bcyl, bcyr;

	rng_t i, j;

	void xchng(int e, int lev = 0) // for previous time levels
	{
          this->xchng(this->mem.psi[e][ this->mem.n[e] - lev], i^halo, j^halo);
	}

	void xchng(arr_2d_t psi, rng_t range_i, rng_t range_j) // for a given array
	{
	  bcxl->fill_halos(psi, range_j);
	  bcxr->fill_halos(psi, range_j);
	  bcyl->fill_halos(psi, range_i);
	  bcyr->fill_halos(psi, range_i);
	}

	// ctor
	solver_2d(mem_t &mem, bc_p &bcxl, bc_p &bcxr, bc_p &bcyl, bc_p &bcyr, const rng_t &i, const rng_t &j) :
	  parent_t(mem),
	  i(i), 
	  j(j),  
	  bcxl(std::move(bcxl)), 
	  bcxr(std::move(bcxr)), 
	  bcyl(std::move(bcyl)),
	  bcyr(std::move(bcyr))
	{}

	public:

	// empty by default
	static void alloc(mem_t &mem, const int nx, const int ny) 
        {
          const rng_t i(0, nx-1), j(0, ny-1);

	  for (int e = 0; e < mem_t::n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem.psi[e].push_back(new arr_2d_t(i^halo, j^halo));

	  mem.C.push_back(new arr_2d_t(i^h   , j^halo));
	  mem.C.push_back(new arr_2d_t(i^halo, j^h   ));
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
