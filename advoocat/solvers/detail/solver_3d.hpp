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
      using namespace arakawa_c;
    
      template<typename real_t, int n_dims, int n_eqs, int n_tlev, int halo>
      class solver_3d : public solver_common<real_t, n_dims, n_eqs, n_tlev, halo>
      {
	using parent_t = solver_common<real_t, n_dims, n_eqs, n_tlev, halo>;

	protected:
      
        typedef std::unique_ptr<bcond::bcond_t<real_t>> bc_p; // TODO: move to parent
	bc_p bcxl, bcxr, bcyl, bcyr, bczl, bczr;

	rng_t i, j, k; // TODO: if stored as idx_t this also could be placed in solver_common

	void xchng(int e, int lev = 0) 
	{
          this->mem->barrier();
	  bcxl->fill_halos(this->mem->psi[e][ this->mem->n[e] - lev ], j^halo, k^halo);
	  bcxr->fill_halos(this->mem->psi[e][ this->mem->n[e] - lev ], j^halo, k^halo);
	  bcyl->fill_halos(this->mem->psi[e][ this->mem->n[e] - lev ], k^halo, i^halo);
	  bcyr->fill_halos(this->mem->psi[e][ this->mem->n[e] - lev ], k^halo, i^halo);
	  bczl->fill_halos(this->mem->psi[e][ this->mem->n[e] - lev ], i^halo, j^halo);
	  bczr->fill_halos(this->mem->psi[e][ this->mem->n[e] - lev ], i^halo, j^halo);
          this->mem->barrier();
	}

	// ctor
	solver_3d(
	  typename parent_t::mem_t *mem,
          bc_p &bcxl,
          bc_p &bcxr,
          bc_p &bcyl,
          bc_p &bcyr,
          bc_p &bczl,
          bc_p &bczr,
	  const rng_t &i,
	  const rng_t &j,
	  const rng_t &k
	) :
	  parent_t(mem),
	  i(i), j(j),  k(k),  
	  bcxl(std::move(bcxl)),
	  bcxr(std::move(bcxr)),
	  bcyl(std::move(bcyl)),
	  bcyr(std::move(bcyr)),
	  bczl(std::move(bczl)),
	  bczr(std::move(bczr))
	{} 

	public:

	static void alloc(
          typename parent_t::mem_t *mem,
          const int nx, const int ny, const int nz
        )   
        {
          const rng_t i(0, nx-1), j(0, ny-1), k(0, nz-1);

	  for (int e = 0; e < n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(new typename parent_t::arr_t(i^halo, j^halo, k^halo)); 

	  mem->C.push_back(new typename parent_t::arr_t(i^h, j^halo, k^halo));
	  mem->C.push_back(new typename parent_t::arr_t(i^halo, j^h, k^halo));
	  mem->C.push_back(new typename parent_t::arr_t(i^halo, j^halo, k^h));
        }  
      };
    };
  }; // namespace solvers
}; // namespace advoocat
