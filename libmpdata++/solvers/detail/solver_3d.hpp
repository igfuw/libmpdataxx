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
      using namespace arakawa_c;
    
      template<typename ct_params_t, int n_tlev, int minhalo>
      class solver<
        ct_params_t, 
        n_tlev, 
        minhalo,
        typename std::enable_if<ct_params_t::n_dims == 3 >::type
      > : public solver_common<ct_params_t, n_tlev, minhalo>
      {
	using parent_t = solver_common<ct_params_t, n_tlev, minhalo>;

	protected:
      
	typename parent_t::bcp_t bcxl, bcxr, bcyl, bcyr, bczl, bczr;

	rng_t i, j, k; // TODO: if stored as idx_t this also could be placed in solver_common
        idx_t<parent_t::n_dims> ijk;

	void xchng(int e, int lev = 0) 
	{
          this->mem->barrier();
	  bcxl->fill_halos_sclr(this->mem->psi[e][ this->n[e] - lev ], j^this->halo, k^this->halo);
	  bcxr->fill_halos_sclr(this->mem->psi[e][ this->n[e] - lev ], j^this->halo, k^this->halo);
	  bcyl->fill_halos_sclr(this->mem->psi[e][ this->n[e] - lev ], k^this->halo, i^this->halo);
	  bcyr->fill_halos_sclr(this->mem->psi[e][ this->n[e] - lev ], k^this->halo, i^this->halo);
	  bczl->fill_halos_sclr(this->mem->psi[e][ this->n[e] - lev ], i^this->halo, j^this->halo);
	  bczr->fill_halos_sclr(this->mem->psi[e][ this->n[e] - lev ], i^this->halo, j^this->halo);
          this->mem->barrier();
	}
        
        void hook_ante_loop(const int nt) // TODO: this nt conflicts in fact with multiple-advance()-call logic!
        {
          parent_t::hook_ante_loop(nt);
	  
          // TODO: make it work with more than one equation    
          this->bcxl->bcinit(this->mem->psi[0][0], this->j, this->k);
          this->bcxr->bcinit(this->mem->psi[0][0], this->j, this->k);
          this->bcyl->bcinit(this->mem->psi[0][0], this->k, this->i);
          this->bcyr->bcinit(this->mem->psi[0][0], this->k, this->i);
          this->bczl->bcinit(this->mem->psi[0][0], this->i, this->j);
          this->bczr->bcinit(this->mem->psi[0][0], this->i, this->j);
          this->bcxl->fill_halos_vctr_alng(this->mem->GC[0], this->j, this->k); // TODO: one xchng?
          this->bcxr->fill_halos_vctr_alng(this->mem->GC[0], this->j, this->k);
          this->bcyl->fill_halos_vctr_alng(this->mem->GC[1], this->k, this->i); // TODO: one xchng?
          this->bcyr->fill_halos_vctr_alng(this->mem->GC[1], this->k, this->i);
          this->bczl->fill_halos_vctr_alng(this->mem->GC[2], this->i, this->j); // TODO: one xchng?
          this->bczr->fill_halos_vctr_alng(this->mem->GC[2], this->i, this->j);
          // sanity check for non-divergence of the initial Courant number field
          // TODO: same in 1D and 3D
          if (blitz::epsilon(typename parent_t::real_t(44)) < max(abs(
	    ( 
              this->mem->GC[0](i-h, j, k) - 
	      this->mem->GC[0](i+h, j, k)
            ) + (
	      this->mem->GC[1](i, j-h, k) - 
	      this->mem->GC[1](i, j+h, k)
            ) + (
	      this->mem->GC[2](i, j, k-h) - 
	      this->mem->GC[2](i, j, k+h)
            )
	  ))) throw std::runtime_error("initial advector field is divergent");
        }

        public:

        struct ctor_args_t
        {   
          typename parent_t::mem_t *mem;
          typename parent_t::bcp_t 
            &bcxl, &bcxr, 
            &bcyl, &bcyr,
            &bczl, &bczr; 
          const rng_t &i, &j, &k; 
        };  

        protected:

	// ctor
	solver(
          ctor_args_t args,
          const typename parent_t::rt_params_t &p
        ) :
	  parent_t(args.mem, p),
	  i(args.i), 
          j(args.j), 
          k(args.k),  
          ijk({args.i, args.j, args.k}),
	  bcxl(std::move(args.bcxl)),
	  bcxr(std::move(args.bcxr)),
	  bcyl(std::move(args.bcyl)),
	  bcyr(std::move(args.bcyr)),
	  bczl(std::move(args.bczl)),
	  bczr(std::move(args.bczr))
	{} 

	public:

	static void alloc(
          typename parent_t::mem_t *mem,
          const typename parent_t::rt_params_t &p
        )   
        {
          // psi
          mem->psi.resize(parent_t::n_eqs);
	  for (int e = 0; e < parent_t::n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t(
                parent_t::rng_sclr(p.span[0]),
                parent_t::rng_sclr(p.span[1]),
                parent_t::rng_sclr(p.span[2])
              ))); 

          // Courant field components (Arakawa-C grid)
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_vctr(p.span[0]),
            parent_t::rng_sclr(p.span[1]),
            parent_t::rng_sclr(p.span[2])
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(p.span[0]),
            parent_t::rng_vctr(p.span[1]),
            parent_t::rng_sclr(p.span[2])
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(p.span[0]),
            parent_t::rng_sclr(p.span[1]),
            parent_t::rng_vctr(p.span[2])
          )));

          // allocate G
          if (formulae::opts::isset(ct_params_t::opts, formulae::opts::nug))
	    mem->G.reset(mem->old(new typename parent_t::arr_t(
                    parent_t::rng_sclr(p.span[0]),
                    parent_t::rng_sclr(p.span[1]),
                    parent_t::rng_sclr(p.span[2])
            )));

	  // allocate Kahan summation temporary vars
	  if (formulae::opts::isset(ct_params_t::opts, formulae::opts::khn))
	    for (int n = 0; n < 3; ++n) 
	      mem->khn_tmp.push_back(mem->old(new typename parent_t::arr_t( 
	        parent_t::rng_sclr(p.span[0]), 
	        parent_t::rng_sclr(p.span[1]),
	        parent_t::rng_sclr(p.span[2])
	      )));
        }  
        
        // helper method to allocate a temporary space composed of vector-component arrays
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem, const std::array<int, 3> &span,
          const char * __file__
        )
        {
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_vctr(span[0]),
            parent_t::rng_sclr(span[1]),
            parent_t::rng_sclr(span[2])
          ))); 
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(span[0]),
            parent_t::rng_vctr(span[1]),
            parent_t::rng_sclr(span[2])
          ))); 
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
             parent_t::rng_sclr(span[0]),
             parent_t::rng_sclr(span[1]),
             parent_t::rng_vctr(span[2])
          ))); 
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem, const std::array<int, 3> &span,
          const char * __file__, const int n_arr
        )   
        {   
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          for (int n = 0; n < n_arr; ++n)
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_sclr(span[0]),
              parent_t::rng_sclr(span[1]),
              parent_t::rng_sclr(span[2])
            )));
        } 
      };
    };
  }; // namespace solvers
}; // namespace libmpdataxx
