/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/solvers/adv/detail/solver_common.hpp>
#include <libmpdata++/arakawa_c.hpp>
#include <libmpdata++/bcond/bcond.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    namespace detail
    {
      using namespace arakawa_c;
    
      template<typename real_t, int n_tlev, formulae::opts::opts_t opts, int minhalo>
      class solver<real_t, 3, n_tlev, opts, minhalo> : public solver_common<real_t, 3, n_tlev, minhalo>
      {
	using parent_t = solver_common<real_t, 3, n_tlev, minhalo>;

	protected:
      
        typedef std::unique_ptr<bcond::bcond_t<real_t>> bc_p; // TODO: move to parent
	bc_p bcxl, bcxr, bcyl, bcyr, bczl, bczr;

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

        public:

        struct ctor_args_t
        {   
          typename parent_t::mem_t *mem;
          bc_p 
            &bcxl, &bcxr, 
            &bcyl, &bcyr,
            &bczl, &bczr; 
          const rng_t &i, &j, &k; 
        };  

        protected:

	// ctor
	solver(
          ctor_args_t args,
          const typename parent_t::params_t &p
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
          const typename parent_t::params_t &p
        )   
        {
	  for (int e = 0; e < p.n_eqs; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t(
                parent_t::rng_sclr(p.span[0]),
                parent_t::rng_sclr(p.span[1]),
                parent_t::rng_sclr(p.span[2])
              ))); 

	  mem->GC.push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_vctr(p.span[0]), parent_t::rng_sclr(p.span[1]), parent_t::rng_sclr(p.span[2])
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(p.span[0]), parent_t::rng_vctr(p.span[1]), parent_t::rng_sclr(p.span[2])
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(p.span[0]), parent_t::rng_sclr(p.span[1]), parent_t::rng_vctr(p.span[2])
          )));

          // TODO: allocate G
        }  
      };
    };
  }; // namespace solvers
}; // namespace libmpdataxx
