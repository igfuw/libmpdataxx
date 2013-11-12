/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/adv/detail/solver_common.hpp>
#include <libmpdata++/bcond/bcond.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      // n_eqs, n_tlev are template parameters as they are needed in both static and non-static contexts
      template<typename real_t, int n_eqs, int n_tlev, int minhalo>
      class solver<real_t, 1, n_eqs, n_tlev, minhalo> : public solver_common<real_t, 1, n_eqs, n_tlev, minhalo>
      {
	using parent_t = solver_common<real_t, 1, n_eqs, n_tlev, minhalo>;

	protected:

        typedef std::unique_ptr<bcond::bcond_t<real_t>> bc_p; // TODO: move to parent
        bc_p bcxl, bcxr;
     
	rng_t i; // TODO: idx_t i do common?
        idx_t<parent_t::n_dims> ijk;

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
	solver(ctor_args_t args) :
	  parent_t(args.mem), 
          bcxl(std::move(args.bcxl)), 
          bcxr(std::move(args.bcxr)),
          i(args.i),
          ijk(args.i)
	{}



        // memory allocation logic using static methods

        private:

	static void alloc_tmp(
	  typename parent_t::mem_t *mem, 
	  const char * __file__, 
	  const int n_arr,
	  const rng_t rng
	)
	{
	  mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>()); 
          for (int n = 0; n < n_arr; ++n)
          {
	    mem->tmp[__file__].back().push_back(
              new typename parent_t::arr_t( rng )
            ); 
          }
	}

        public:

	static void alloc(typename parent_t::mem_t *mem, const int nx)   
        {
	  for (int e = 0; e < n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(new typename parent_t::arr_t(parent_t::rng_sclr(nx)));
    
	  mem->GC.push_back(new typename parent_t::arr_t(parent_t::rng_vctr(nx))); 

          // TODO: allocate G
          assert(mem->G.numElements() == 0);
        } 

        protected:

        // helper method to allocate a vector-component temporary array
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem, const int nx, 
          const char * __file__
        )
        {
          alloc_tmp(mem, __file__, 1, parent_t::rng_vctr(nx)); // always one-component vectors
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem, const int nx, 
          const char * __file__, const int n_arr
        )
        {
          alloc_tmp(mem, __file__, n_arr, parent_t::rng_sclr(nx)); 
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
