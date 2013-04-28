/** @file 
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <advoocat/solvers/detail/solver_common.hpp>
#include <advoocat/arakawa_c.hpp>
#include <advoocat/bcond/bcond.hpp>

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      using namespace advoocat::arakawa_c;

      template<typename real_t, int n_dims, int n_eqs, int n_tlev, int halo>
      class solver_2d : public solver_common<real_t, n_dims, n_eqs, n_tlev, halo>
      {
	using parent_t = solver_common<real_t, n_dims, n_eqs, n_tlev, halo>;

	protected:
      
        typedef std::unique_ptr<bcond::bcond_t<real_t>> bc_p; // TODO: move to parent
	bc_p bcxl, bcxr, bcyl, bcyr;

	rng_t i, j;
        idx_t<n_dims> ijk;

	void xchng(int e, int lev = 0) // for previous time levels
	{
          this->xchng(this->mem->psi[e][ this->n[e] - lev], i^halo, j^halo);
	}

	void xchng(typename parent_t::arr_t psi, rng_t range_i, rng_t range_j) // for a given array
	{
          this->mem->barrier();
	  bcxl->fill_halos(psi, range_j);
	  bcxr->fill_halos(psi, range_j);
	  bcyl->fill_halos(psi, range_i);
	  bcyr->fill_halos(psi, range_i);
          this->mem->barrier();
	}

        public:
 
        struct ctor_args_t
        {   
          typename parent_t::mem_t *mem;
          bc_p &bcxl, &bcxr, &bcyl, &bcyr; 
          const rng_t &i, &j; 
        };  

        protected:

	// ctor
	solver_2d(ctor_args_t args) :
	  parent_t(args.mem),
	  i(args.i), 
	  j(args.j),  
          ijk({args.i, args.j}),
	  bcxl(std::move(args.bcxl)), 
	  bcxr(std::move(args.bcxr)), 
	  bcyl(std::move(args.bcyl)),
	  bcyr(std::move(args.bcyr))
	{}

	public:

	static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny) 
        {
          const rng_t i(0, nx-1), j(0, ny-1);

	  for (int e = 0; e < n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(new typename parent_t::arr_t(i^halo, j^halo));

	  mem->C.push_back(new typename parent_t::arr_t(i^h   , j^halo));
	  mem->C.push_back(new typename parent_t::arr_t(i^halo, j^h   ));
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace advoocat
