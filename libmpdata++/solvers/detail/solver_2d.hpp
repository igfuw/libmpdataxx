/** @file 
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/solvers/detail/solver_common.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      using namespace libmpdataxx::arakawa_c;

      template<typename ct_params_t, int n_tlev, int minhalo>
      class solver<
        ct_params_t, 
        n_tlev,  
        minhalo,
        typename std::enable_if<ct_params_t::n_dims == 2 >::type
      > : public solver_common<ct_params_t, n_tlev, minhalo>
      {
	using parent_t = solver_common<ct_params_t, n_tlev, minhalo>;

	protected:
      
	typename parent_t::bcp_t bcxl, bcxr, bcyl, bcyr;

	const rng_t i, j; // TODO: to be removed

	void xchng_sclr(typename parent_t::arr_t arr, rng_t range_i, rng_t range_j, const bool deriv = false) // for a given array
	{
          this->mem->barrier();
	  bcxl->fill_halos_sclr(arr, range_j, deriv);
	  bcxr->fill_halos_sclr(arr, range_j, deriv);
	  bcyl->fill_halos_sclr(arr, range_i, deriv);
	  bcyr->fill_halos_sclr(arr, range_i, deriv);
          this->mem->barrier();
	}

	void xchng(int e) 
	{
          this->xchng_sclr(this->mem->psi[e][ this->n[e]], i^this->halo, j^this->halo);
	}

        virtual void xchng_vctr_alng(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          this->mem->barrier();
          bcxl->fill_halos_vctr_alng(arrvec, j);
          bcxr->fill_halos_vctr_alng(arrvec, j);
          bcyl->fill_halos_vctr_alng(arrvec, i);
          bcyr->fill_halos_vctr_alng(arrvec, i);
          // TODO: open bc nust be last!!!
          this->mem->barrier();
        }

        // TODO: ref in argument...
        void hook_ante_loop(const int nt) // TODO: this nt conflicts in fact with multiple-advance()-call logic!
        {
          parent_t::hook_ante_loop(nt);

          xchng_vctr_alng(this->mem->GC);

          // sanity check for non-divergence of the initial Courant number field
          // (including compatibility with the initial condition)
          // TODO: same in 1D
          {
            typename ct_params_t::real_t max_abs_div = max(abs(
	      ( 
		this->mem->GC[0](i-h, j  ) - 
		this->mem->GC[0](i+h, j  )
	      ) + (
		this->mem->GC[1](i,   j-h) - 
		this->mem->GC[1](i,   j+h)
	      )
	    ));
	    if (max_abs_div > blitz::epsilon(typename parent_t::real_t(44))) 
	      throw std::runtime_error("initial advector field is divergent");
          }
        }

        public:
 
        struct ctor_args_t
        {   
          typename parent_t::mem_t *mem;
          typename parent_t::bcp_t &bcxl, &bcxr, &bcyl, &bcyr; 
          const rng_t &i, &j; 
        };  
        
        struct rt_params_t : parent_t::rt_params_t
        {
          typename parent_t::real_t di = 0, dj = 0;
        };

        protected:

	// ctor
	solver(
          ctor_args_t args,
          const rt_params_t &p
        ) :
	  parent_t(
            args.mem, 
            p, 
            idx_t<parent_t::n_dims>({args.i, args.j})
          ),
	  i(args.i), 
	  j(args.j),  
	  bcxl(std::move(args.bcxl)), 
	  bcxr(std::move(args.bcxr)), 
	  bcyl(std::move(args.bcyl)),
	  bcyr(std::move(args.bcyr))
	{
          this->di = p.di;
          this->dj = p.dj;
        }

        // memory allocation logic using static methods

	public:

	static void alloc(typename parent_t::mem_t *mem, const rt_params_t &p) 
        {
          // psi 
          mem->psi.resize(parent_t::n_eqns);
	  for (int e = 0; e < parent_t::n_eqns; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t( 
                parent_t::rng_sclr(p.grid_size[0]), 
                parent_t::rng_sclr(p.grid_size[1])
              )));

          // Courant field components (Arakawa-C grid)
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_vctr(p.grid_size[0]), 
            parent_t::rng_sclr(p.grid_size[1]) 
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_sclr(p.grid_size[0]), 
            parent_t::rng_vctr(p.grid_size[1]) 
          )));
 
          // allocate G
          if (opts::isset(ct_params_t::opts, opts::nug))
	    mem->G.reset(mem->old(new typename parent_t::arr_t(
                    parent_t::rng_sclr(p.grid_size[0]),
                    parent_t::rng_sclr(p.grid_size[1])
            )));

          // allocate Kahan summation temporary vars
          if (!opts::isset(ct_params_t::opts, opts::nkh))
	    for (int n = 0; n < 3; ++n) 
	      mem->khn_tmp.push_back(mem->old(new typename parent_t::arr_t( 
                parent_t::rng_sclr(p.grid_size[0]), 
                parent_t::rng_sclr(p.grid_size[1])
              )));
        }

        protected:

        // helper method to allocate a temporary space composed of vector-component arrays
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem, const std::array<int, 2> &grid_size,
          const char * __file__
        )
        {
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_vctr(grid_size[0]), 
            parent_t::rng_sclr(grid_size[1]) 
          ))); 
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_sclr(grid_size[0]), 
            parent_t::rng_vctr(grid_size[1]) 
          ))); 
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem, const std::array<int, 2> &grid_size,
          const char * __file__, const int n_arr
        )   
        {   
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          for (int n = 0; n < n_arr; ++n)
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_sclr(grid_size[0]),
              parent_t::rng_sclr(grid_size[1])
            )));
        } 
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
