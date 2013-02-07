/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "solver_common.hpp"
#include "../arakawa_c.hpp"

namespace advoocat
{
  namespace solvers
  {
    namespace detail
    {
      using namespace arakawa_c;
    
      template<class bcx_t, class bcy_t, class bcz_t, class mem_t, int n_tlev>
      class solver_3d : public solver_common<mem_t, n_tlev>
      {
	using parent_t = solver_common<mem_t, n_tlev>;
	using arr_3d_t = typename mem_t::arr_t;

	protected:
      
	bcx_t bcx;
	bcy_t bcy;
	bcz_t bcz;

	// <TODO> this could in principle be placed in solver_common
	int halo;
	// </TODO>

	rng_t i, j, k; // TODO: if stored as idx_t this also could be placed in solver_common

	void xchng(int e, int lev = 0) 
	{
	  bcx.fill_halos(this->mem.psi[e][ this->mem.n[e] - lev ], j^halo, k^halo);
	  bcy.fill_halos(this->mem.psi[e][ this->mem.n[e] - lev ], k^halo, i^halo);
	  bcz.fill_halos(this->mem.psi[e][ this->mem.n[e] - lev ], i^halo, j^halo);
	}

	// ctor
	solver_3d(
	  mem_t &mem,
	  const rng_t &i,
	  const rng_t &j,
	  const rng_t &k,
	  int halo
	) :
	  parent_t(mem),
	  halo(halo),
	  i(i), j(j),  k(k),  
	  bcx(rng_t(i), halo), 
	  bcy(rng_t(j), halo),
	  bcz(rng_t(k), halo)
	{} 

	public:

	// empty by default
	static void alloc(mem_t &mem, const int nx, const int ny, const int nz)   
        {
          const rng_t i(0, nx-1), j(0, ny-1), k(0, nz-1);
          const int hlo = 1; // TODO!!

	  for (int e = 0; e < mem_t::n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem.psi[e].push_back(new arr_3d_t(i^hlo, j^hlo, k^hlo)); 

	  mem.C.push_back(new arr_3d_t(i^h, j^hlo, k^hlo));
	  mem.C.push_back(new arr_3d_t(i^hlo, j^h, k^hlo));
	  mem.C.push_back(new arr_3d_t(i^hlo, j^hlo, k^h));
        }  
      };
    };
  }; // namespace solvers
}; // namespace advoocat
