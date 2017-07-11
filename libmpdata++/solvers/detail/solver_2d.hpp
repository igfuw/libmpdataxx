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
      
	const rng_t i, j; // TODO: to be removed
        typename parent_t::arr_t &courant_field; // TODO: should be in solver common but cannot be allocated there ?

	virtual void xchng_sclr(typename parent_t::arr_t &arr,
                        const idx_t<2> &range_ijk,
                        const int ext = 0,
                        const bool deriv = false
        ) final // for a given array
	{
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_sclr(arr, range_ijk[1]^ext, deriv);
	  for (auto &bc : this->bcs[1]) bc->fill_halos_sclr(arr, range_ijk[0]^ext, deriv);
          this->mem->barrier();
	}

	void xchng(int e) final
	{
          this->xchng_sclr(this->mem->psi[e][ this->n[e]], this->ijk, this->halo);
	}

        void xchng_vctr_alng(const arrvec_t<typename parent_t::arr_t> &arrvec, const typename parent_t::real_t flux = 0) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_alng(arrvec, j, flux);
          for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_alng(arrvec, i, flux);
          // TODO: open bc nust be last!!!
          this->mem->barrier();
        }

        virtual void xchng_vctr_nrml(
          const arrvec_t<typename parent_t::arr_t> &arrvec, 
          const idx_t<2> &range_ijk,
          const int ext = 0
        ) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_nrml(arrvec[0], range_ijk[0]^ext^h);
          for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_nrml(arrvec[1], range_ijk[1]^ext^h);
          this->mem->barrier();
        }

        virtual void xchng_pres(
          const typename parent_t::arr_t &arr,
          const idx_t<2> &range_ijk
        ) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_pres(arr, range_ijk[1]);
          for (auto &bc : this->bcs[1]) bc->fill_halos_pres(arr, range_ijk[0]);
          this->mem->barrier();
        }

        virtual void set_edges(
          const arrvec_t<typename parent_t::arr_t> &av,
          const idx_t<2> &range_ijk,
          const int &sign
        ) final
        {
          for (auto &bc : this->bcs[0]) bc->set_edge_pres(av[0], range_ijk[1], sign);
          for (auto &bc : this->bcs[1]) bc->set_edge_pres(av[1], range_ijk[0], sign);
          this->mem->barrier();
        }

        virtual void save_edges(
          const arrvec_t<typename parent_t::arr_t> &av,
          const idx_t<2> &range_ijk
        ) final
        {
          for (auto &bc : this->bcs[0]) bc->save_edge_vel(av[0], range_ijk[1]);
          for (auto &bc : this->bcs[1]) bc->save_edge_vel(av[1], range_ijk[0]);
          this->mem->barrier();
        }

        // TODO: ref in argument...
        void hook_ante_loop(const typename parent_t::advance_arg_t nt) // TODO: this nt conflicts in fact with multiple-advance()-call logic!
        {
          parent_t::hook_ante_loop(nt);

          // sanity check for non-divergence of the initial Courant number field
          // (including compatibility with the initial condition)
          // TODO: same in 1D
          if (!opts::isset(ct_params_t::opts, opts::dfl))
          {
            typename ct_params_t::real_t max_abs_div = max(abs(
              (
                ( 
                  this->mem->GC[0](i-h, j  ) - 
                  this->mem->GC[0](i+h, j  )
                ) + (
                  this->mem->GC[1](i,   j-h) - 
                  this->mem->GC[1](i,   j+h)
                )
              ) / formulae::G<ct_params_t::opts, 0>(*this->mem->G, i, j)
	    ));

	    if (max_abs_div > this->max_abs_div_eps) 
	      throw std::runtime_error("initial advector field is divergent");
          }
        }

        typename parent_t::real_t courant_number(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          courant_field(this->ijk) = 0.5 * (
                                             abs(arrvec[0](i+h, j) + arrvec[0](i-h, j))
                                           + abs(arrvec[1](i, j+h) + arrvec[1](i, j-h))
                                           ) / formulae::G<ct_params_t::opts, 0>(*this->mem->G, i, j);
          return this->mem->max(this->rank, courant_field(this->ijk));
        }
        
        void scale_gc(const typename parent_t::real_t time,
                      const typename parent_t::real_t cur_dt,
                      const typename parent_t::real_t old_dt) final
        {
          this->mem->GC[0](rng_t(i.first(), i.last()-1)^h, j) *= cur_dt / old_dt;
          this->mem->GC[1](i, rng_t(j.first(), j.last()-1)^h) *= cur_dt / old_dt;
          this->xchng_vctr_alng(this->mem->GC);
          auto ex = this->halo - 1;
          this->xchng_vctr_nrml(this->mem->GC, this->ijk, ex);
        }

        public:
 
        struct ctor_args_t
        {   
          // <TODO> these should be common for 1D,2D,3D
          int rank;
          typename parent_t::mem_t *mem;
          // </TODO>
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
            args.rank,
            args.mem, 
            p, 
            idx_t<parent_t::n_dims>({args.i, args.j})
          ),
	  i(args.i), 
	  j(args.j),
          courant_field(args.mem->tmp[__FILE__][0][0])
	{
          this->di = p.di;
          this->dj = p.dj;
          this->dijk = {p.di, p.dj};
	  this->set_bcs(0, args.bcxl, args.bcxr); 
	  this->set_bcs(1, args.bcyl, args.bcyr);
        }

        // memory allocation logic using static methods

	public:

	static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
          // psi 
          mem->psi.resize(parent_t::n_eqns);
	  for (int e = 0; e < parent_t::n_eqns; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t( 
                parent_t::rng_sclr(mem->grid_size[0]), 
                parent_t::rng_sclr(mem->grid_size[1])
              )));

          // Courant field components (Arakawa-C grid)
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_vctr(mem->grid_size[0]), 
            parent_t::rng_sclr(mem->grid_size[1]) 
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_sclr(mem->grid_size[0]), 
            parent_t::rng_vctr(mem->grid_size[1]) 
          )));

          // fully third-order accurate mpdata needs also time derivatives of
          // the Courant field
          if (opts::isset(ct_params_t::opts, opts::div_3rd))
          {
            // TODO: why for (auto f : {mem->ndt_GC, mem->ndtt_GC}) doesn't work ?
            mem->ndt_GC.push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_vctr(mem->grid_size[0]), 
              parent_t::rng_sclr(mem->grid_size[1]) 
            )));
            mem->ndt_GC.push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_sclr(mem->grid_size[0]), 
              parent_t::rng_vctr(mem->grid_size[1]) 
            )));
            mem->ndtt_GC.push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_vctr(mem->grid_size[0]), 
              parent_t::rng_sclr(mem->grid_size[1]) 
            )));
            mem->ndtt_GC.push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_sclr(mem->grid_size[0]), 
              parent_t::rng_vctr(mem->grid_size[1]) 
            )));
          }
 
          // allocate G
          if (opts::isset(ct_params_t::opts, opts::nug))
	    mem->G.reset(mem->old(new typename parent_t::arr_t(
                    parent_t::rng_sclr(mem->grid_size[0]),
                    parent_t::rng_sclr(mem->grid_size[1])
            )));

          // allocate Kahan summation temporary vars
          if (opts::isset(ct_params_t::opts, opts::khn))
	    for (int n = 0; n < 3; ++n) 
	      mem->khn_tmp.push_back(mem->old(new typename parent_t::arr_t( 
                parent_t::rng_sclr(mem->grid_size[0]), 
                parent_t::rng_sclr(mem->grid_size[1])
              )));
          // courant field
          alloc_tmp_sclr(mem, __FILE__, 1);
        }

        protected:

        // helper method to allocate a temporary space composed of vector-component arrays
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem,
          const char * __file__
        )
        {
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_vctr(mem->grid_size[0]), 
            parent_t::rng_sclr(mem->grid_size[1]) 
          ))); 
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_sclr(mem->grid_size[0]), 
            parent_t::rng_vctr(mem->grid_size[1]) 
          ))); 
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem, 
          const char * __file__, const int n_arr,
          std::string name = "",
          bool srfc = false
        )   
        {   
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());

          if (!name.empty()) mem->avail_tmp[name] = std::make_pair(__file__, mem->tmp[__file__].size() - 1);

          for (int n = 0; n < n_arr; ++n)
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_sclr(mem->grid_size[0]),
              srfc ? rng_t(0, 0) : parent_t::rng_sclr(mem->grid_size[1])
            )));
        } 
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
