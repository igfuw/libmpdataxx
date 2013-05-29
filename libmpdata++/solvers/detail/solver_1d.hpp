/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/detail/solver_common.hpp>
#include <libmpdata++/bcond/bcond.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      using namespace libmpdataxx::arakawa_c;

      template<typename real_t, int n_dims, int n_eqs, int n_tlev, int halo>
      class solver_1d : public solver_common<real_t, n_dims, n_eqs, n_tlev, halo>
      {
	using parent_t = solver_common<real_t, n_dims, n_eqs, n_tlev, halo>;

	protected:

        typedef std::unique_ptr<bcond::bcond_t<real_t>> bc_p; // TODO: move to parent
        bc_p bcxl, bcxr;
     
	rng_t i; // TODO: idx_t i do common?
        idx_t<n_dims> ijk;

	void xchng(int e, int lev = 0) 
	{
          this->mem->barrier();
	  bcxl->fill_halos_sclr( this->mem->psi[e][ this->n[e] - lev ] );
	  bcxr->fill_halos_sclr( this->mem->psi[e][ this->n[e] - lev ] );
          this->mem->barrier();
	}

        public:
 
        struct ctor_args_t
        {
          typename parent_t::mem_t *mem;
          bc_p &bcxl, &bcxr; 
          const rng_t &i;
        };

        protected:

	// ctor
	solver_1d(ctor_args_t args) :
	  parent_t(args.mem), 
          bcxl(std::move(args.bcxl)), 
          bcxr(std::move(args.bcxr)),
          i(args.i),
          ijk(args.i)
	{}

	public:

	static void alloc(typename parent_t::mem_t *mem, const int nx)   
        {
// TODO: N_DEBUG code to assure all parents called (same in 2D, 3D)
          const rng_t i(0, nx-1);
          
	  for (int e = 0; e < n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(new typename parent_t::arr_t(i^halo));
    
	  mem->C.push_back(new typename parent_t::arr_t(i^h^(halo-1)));
        } 
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
