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
      template<typename ct_params_t, int n_tlev, int minhalo>
      class solver<
        ct_params_t, 
        n_tlev, 
        minhalo,
        typename std::enable_if<ct_params_t::n_dims == 1 >::type
      > : public solver_common<ct_params_t, n_tlev, minhalo>
      {
	using parent_t = solver_common<ct_params_t, n_tlev, minhalo>;

	protected:

        typedef std::unique_ptr<bcond::bcond_t<typename parent_t::real_t>> bc_p; // TODO: move to parent
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
	solver(
          ctor_args_t args,
          const typename parent_t::rt_params_t &p
        ) :
	  parent_t(args.mem, p), 
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
              mem->old(new typename parent_t::arr_t( rng ))
            ); 
          }
	}

        public:

	static void alloc(typename parent_t::mem_t *mem, const typename parent_t::rt_params_t &p)   
        {
          mem->psi.resize(parent_t::n_eqs);
	  for (int e = 0; e < parent_t::n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_sclr(p.span[0]))));
    
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_vctr(p.span[0])))); 

          if (formulae::opts::isset(ct_params_t::opts, formulae::opts::nug))
	    mem->G.reset(mem->old(new typename parent_t::arr_t(parent_t::rng_sclr(p.span[0]))));
        } 

        protected:

        // helper method to allocate a vector-component temporary array
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem, const std::array<int, 1> &span, 
          const char * __file__
        )
        {
          alloc_tmp(mem, __file__, 1, parent_t::rng_vctr(span[0])); // always one-component vectors
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem, const std::array<int, 1> &span, 
          const char * __file__, const int n_arr
        )
        {
          alloc_tmp(mem, __file__, n_arr, parent_t::rng_sclr(span[0])); 
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
