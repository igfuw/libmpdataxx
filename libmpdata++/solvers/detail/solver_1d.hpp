/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/opts.hpp>
#include <libmpdata++/solvers/detail/solver_common.hpp>

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
        
        public:

        using typename parent_t::real_t;

	protected:

        typename parent_t::bcp_t bcxl, bcxr;
     
	const rng_t i; //TODO: to be removed

	void xchng(int e) 
	{
          this->mem->barrier(); // TODO: implement using the xchng below
	  bcxl->fill_halos_sclr( this->mem->psi[e][ this->n[e]] );
	  bcxr->fill_halos_sclr( this->mem->psi[e][ this->n[e]] );
          this->mem->barrier();
	}

        virtual void xchng_sclr(typename parent_t::arr_t &arr, const bool deriv = false) final // for a given array
        {
          this->mem->barrier();
          bcxl->fill_halos_sclr(arr, deriv);
          bcxr->fill_halos_sclr(arr, deriv);
          this->mem->barrier();
        }

        virtual void xchng_vctr_alng(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          this->mem->barrier();
          bcxl->fill_halos_vctr_alng(arrvec); 
          bcxr->fill_halos_vctr_alng(arrvec);
          this->mem->barrier();
        }

        void hook_ante_loop(const real_t tshift) // TODO: this tshift conflicts in fact with multiple-advance()-call logic!
        {
          parent_t::hook_ante_loop(tshift);
	  
          // filling halo in velocity field
          xchng_vctr_alng(this->mem->GC);
        }

        public:
 
        struct ctor_args_t
        {
          // <TODO> these should be common for 1D,2D,3D
          int rank;
          typename parent_t::mem_t *mem;
          // </TODO>
          typename parent_t::bcp_t &bcxl, &bcxr; 
          const rng_t &i;
        };

        struct rt_params_t : parent_t::rt_params_t
        {
          real_t di = 0;
        };

        protected:

	// ctor
	solver(
          ctor_args_t args,
          const rt_params_t &p
        ) :
	  parent_t(
            args.rank,
            args.mem, 
            p, 
            idx_t<parent_t::n_dims>(args.i)
          ), 
          bcxl(std::move(args.bcxl)), 
          bcxr(std::move(args.bcxr)),
          i(args.i)
	{
          this->di = p.di;
        }

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

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p)   
        {
          mem->psi.resize(parent_t::n_eqns);
	  for (int e = 0; e < parent_t::n_eqns; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_sclr(p.grid_size[0]))));
    
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_vctr(p.grid_size[0])))); 

          if (opts::isset(ct_params_t::opts, opts::nug))
	    mem->G.reset(mem->old(new typename parent_t::arr_t(parent_t::rng_sclr(p.grid_size[0]))));

          // allocate Kahan summation temporary vars
          if (opts::isset(ct_params_t::opts, opts::khn))
            for (int n = 0; n < 3; ++n) 
              mem->khn_tmp.push_back(mem->old(new typename parent_t::arr_t( 
                parent_t::rng_sclr(p.grid_size[0])
              )));
        } 

        protected:

        // helper method to allocate a vector-component temporary array
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem, const std::array<int, 1> &grid_size, 
          const char * __file__
        )
        {
          alloc_tmp(mem, __file__, 1, parent_t::rng_vctr(grid_size[0])); // always one-component vectors
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem, const std::array<int, 1> &grid_size, 
          const char * __file__, const int n_arr
        )
        {
          alloc_tmp(mem, __file__, n_arr, parent_t::rng_sclr(grid_size[0])); 
        }
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
