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

	const rng_t i, j, k; // TODO: we have ijk in solver_common - could it be removed?

        // generic field used for various statistics (currently Courant number and divergence)
        typename parent_t::arr_t &stat_field; // TODO:/: should be in solver common but cannot be allocated there ?

	virtual void xchng_sclr(typename parent_t::arr_t &arr,
	               const idx_t<3> &range_ijk,
                       const int ext = 0,
                       const bool deriv = false
        ) final // for a given array
	{
          const auto range_ijk_0__ext = this->extend_range(range_ijk[0], ext);
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_sclr(arr, range_ijk[1]^ext, range_ijk[2]^ext, deriv);
	  for (auto &bc : this->bcs[1]) bc->fill_halos_sclr(arr, range_ijk[2]^ext, range_ijk_0__ext, deriv);
	  for (auto &bc : this->bcs[2]) bc->fill_halos_sclr(arr, range_ijk_0__ext, range_ijk[1]^ext, deriv);
          this->mem->barrier();
	}
	void xchng(int e) final
	{
	  this->xchng_sclr(this->mem->psi[e][ this->n[e]], this->ijk, this->halo);
	}

        void xchng_vctr_alng(arrvec_t<typename parent_t::arr_t> &arrvec, const bool ad = false, const bool cyclic = false) final
        {
          this->mem->barrier();
          if (!cyclic)
          {
            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_alng(arrvec, j, k, ad); 
            for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_alng(arrvec, k, i, ad); 
            for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_alng(arrvec, i, j, ad);
          }
          else
          {
            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_alng_cyclic(arrvec, j, k, ad); 
            for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_alng_cyclic(arrvec, k, i, ad); 
            for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_alng_cyclic(arrvec, i, j, ad);
          }
          this->mem->barrier();
        }
        
        virtual void xchng_flux(arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_flux(arrvec, j, k);
          for (auto &bc : this->bcs[1]) bc->fill_halos_flux(arrvec, k, i);
          for (auto &bc : this->bcs[2]) bc->fill_halos_flux(arrvec, i, j);
        }
        
        virtual void xchng_sgs_div(
	  typename parent_t::arr_t &arr,
	  const idx_t<3> &range_ijk
        ) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_sgs_div(arr, range_ijk[1], range_ijk[2]^h);
          for (auto &bc : this->bcs[1]) bc->fill_halos_sgs_div(arr, range_ijk[2]^h, range_ijk[0]);
          for (auto &bc : this->bcs[2]) bc->fill_halos_sgs_div(arr, range_ijk[0], range_ijk[1]);
          this->mem->barrier();
        }
	
        virtual void xchng_sgs_vctr(arrvec_t<typename parent_t::arr_t> &av,
                                    const typename parent_t::arr_t &b, 
	                            const idx_t<3> &range_ijk
        ) final
	{
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_sgs_vctr(av, b, range_ijk[1], range_ijk[2]);
	  for (auto &bc : this->bcs[1]) bc->fill_halos_sgs_vctr(av, b, range_ijk[2], range_ijk[0]);
	  for (auto &bc : this->bcs[2]) bc->fill_halos_sgs_vctr(av, b, range_ijk[0], range_ijk[1]);
          this->mem->barrier();
	}
        
        virtual void xchng_sgs_tnsr_diag(arrvec_t<typename parent_t::arr_t> &av,
                                         const typename parent_t::arr_t &w,
                                         const typename parent_t::arr_t &vip_div,
	                                 const idx_t<3> &range_ijk
        ) final
	{
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_sgs_tnsr(av, w, vip_div, range_ijk[1], range_ijk[2], this->dijk[0]);
	  for (auto &bc : this->bcs[1]) bc->fill_halos_sgs_tnsr(av, w, vip_div, range_ijk[2], range_ijk[0], this->dijk[1]);
	  for (auto &bc : this->bcs[2]) bc->fill_halos_sgs_tnsr(av, w, vip_div, range_ijk[0], range_ijk[1], this->dijk[2]);
          this->mem->barrier();
	}
        
        virtual void xchng_sgs_tnsr_offdiag(arrvec_t<typename parent_t::arr_t> &av,
                                            const arrvec_t<typename parent_t::arr_t> &bv, 
	                                    const idx_t<3> &range_ijk,
	                                    const std::array<rng_t, 3> &range_ijkm
        ) final
	{
          // off-diagonal components of stress tensor are treated the same as a vector
          this->mem->barrier();
          for (auto &bc : this->bcs[0])
          {
            bc->fill_halos_sgs_vctr(av, bv[0], range_ijkm[1], range_ijk[2]^1, 3);
            bc->fill_halos_sgs_vctr(av, bv[1], range_ijk[1]^1, range_ijkm[2], 4);
          }

	  for (auto &bc : this->bcs[1])
          {
            bc->fill_halos_sgs_vctr(av, bv[0], range_ijk[2]^1, range_ijkm[0], 2);
            bc->fill_halos_sgs_vctr(av, bv[1], range_ijkm[2], range_ijk[0]^1, 4);
          }

	  for (auto &bc : this->bcs[2])
          {
            bc->fill_halos_sgs_vctr(av, bv[0], range_ijkm[0], range_ijk[1]^1, 2);
            bc->fill_halos_sgs_vctr(av, bv[1], range_ijk[0]^1, range_ijkm[1], 3);
          }
          this->mem->barrier();
	}

        virtual void xchng_vctr_nrml(
          arrvec_t<typename parent_t::arr_t> &arrvec,
          const idx_t<3> &range_ijk,
          const int ext = 0,
          const bool cyclic = false
        ) final
        {
          this->mem->barrier();
          const auto range_ijk_0__ext_h = this->extend_range(range_ijk[0], ext, h);
          const auto range_ijk_0__ext_1 = this->extend_range(range_ijk[0], ext, 1);
          if (!cyclic)
          {
            for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_nrml(arrvec[0], range_ijk[2]^ext^1, range_ijk[0]^ext^h);
  
            // without this barrier, there is a race condition when some threads handle subdomains
            // with one gridpoint width, the problem manifests itself, for example, in pbl test
            // TODO: figure out the exact cause and try to avoid this barrier, what about 2D,
            //       what about atypical boundary condition choices -- rigid/cyclic/rigid etc
            if (parent_t::div3_mpdata)
            {
              this->mem->barrier();
            }

            for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_nrml(arrvec[0], range_ijk_0__ext_h, range_ijk[1]^ext^1);

            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_nrml(arrvec[1], range_ijk[1]^ext^h, range_ijk[2]^ext^1);
            for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_nrml(arrvec[1], range_ijk_0__ext_1, range_ijk[1]^ext^h);
       
            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_nrml(arrvec[2], range_ijk[1]^ext^1, range_ijk[2]^ext^h);
            for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_nrml(arrvec[2], range_ijk[2]^ext^h, range_ijk_0__ext_1);
          }
          else
          {
            for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_nrml_cyclic(arrvec[0], range_ijk[2]^ext^1, range_ijk_0__ext_h);
            for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_nrml_cyclic(arrvec[0], range_ijk_0__ext_h, range_ijk[1]^ext^1);

            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_nrml_cyclic(arrvec[1], range_ijk[1]^ext^h, range_ijk[2]^ext^1);
            for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_nrml_cyclic(arrvec[1], range_ijk_0__ext_1, range_ijk[1]^ext^h);
       
            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_nrml_cyclic(arrvec[2], range_ijk[1]^ext^1, range_ijk[2]^ext^h);
            for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_nrml_cyclic(arrvec[2], range_ijk[2]^ext^h, range_ijk_0__ext_1);
          }
          this->mem->barrier();
        }

        virtual void xchng_pres(
	  typename parent_t::arr_t &arr,
	  const idx_t<3> &range_ijk,
          const int ext = 0
        ) final
        {
          const auto range_ijk_0__ext = this->extend_range(range_ijk[0], ext);
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_pres(arr, range_ijk[1]^ext, range_ijk[2]^ext);
          for (auto &bc : this->bcs[1]) bc->fill_halos_pres(arr, range_ijk[2]^ext, range_ijk_0__ext);
          for (auto &bc : this->bcs[2]) bc->fill_halos_pres(arr, range_ijk_0__ext, range_ijk[1]^ext);
          this->mem->barrier();
        }

        virtual void set_edges(
          arrvec_t<typename parent_t::arr_t> &av,
	  const idx_t<3> &range_ijk,
          const int &sign
        ) final
        {
          for (auto &bc : this->bcs[0]) bc->set_edge_pres(av[0], range_ijk[1], range_ijk[2], sign);
          for (auto &bc : this->bcs[1]) bc->set_edge_pres(av[1], range_ijk[2], range_ijk[0], sign);
          for (auto &bc : this->bcs[2]) bc->set_edge_pres(av[2], range_ijk[0], range_ijk[1], sign);
          this->mem->barrier();
        }
        
        virtual void save_edges(
          const arrvec_t<typename parent_t::arr_t> &av,
	  const idx_t<3> &range_ijk
        ) final
        {
          for (auto &bc : this->bcs[0]) bc->save_edge_vel(av[0], range_ijk[1], range_ijk[2]);
          for (auto &bc : this->bcs[1]) bc->save_edge_vel(av[1], range_ijk[2], range_ijk[0]);
          for (auto &bc : this->bcs[2]) bc->save_edge_vel(av[2], range_ijk[0], range_ijk[1]);
          this->mem->barrier();
        }
        
        void hook_ante_loop(const typename parent_t::advance_arg_t nt) // TODO: this nt conflicts in fact with multiple-advance()-call logic!
        {
          parent_t::hook_ante_loop(nt);
 
          // sanity check for non-divergence of the initial Courant number field
          // TODO: same in 1D
          if (!opts::isset(ct_params_t::opts, opts::dfl))
          {
            typename ct_params_t::real_t max_abs_div = max_abs_vctr_div(this->mem->GC);

	    if (max_abs_div > this->max_abs_div_eps) 
	      throw std::runtime_error("initial advector field is divergent");
          }
        }

        typename parent_t::real_t courant_number(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          stat_field(this->ijk) =  0.5 * (
                                            abs(arrvec[0](i+h, j, k) + arrvec[0](i-h, j, k))
                                          + abs(arrvec[1](i, j+h, k) + arrvec[1](i, j-h, k))
                                          + abs(arrvec[2](i, j, k+h) + arrvec[2](i, j, k-h))
                                         ) / formulae::G<ct_params_t::opts, 0>(*this->mem->G, i, j, k);
          return this->mem->max(this->rank, stat_field(this->ijk));
        }
        
        typename parent_t::real_t max_abs_vctr_div(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          stat_field(this->ijk) =  abs(
                                         (arrvec[0](i+h, j, k) - arrvec[0](i-h, j, k))
                                       + (arrvec[1](i, j+h, k) - arrvec[1](i, j-h, k))
                                       + (arrvec[2](i, j, k+h) - arrvec[2](i, j, k-h))
                                      ) / formulae::G<ct_params_t::opts, 0>(*this->mem->G, i, j, k);

          return this->mem->max(this->rank, stat_field(this->ijk));
        }

        void scale_gc(const typename parent_t::real_t time,
                      const typename parent_t::real_t cur_dt,
                      const typename parent_t::real_t old_dt) final
        {
          this->mem->GC[0](rng_t(i.first(), i.last()-1)^h, j, k) *= cur_dt / old_dt;
          this->mem->GC[1](i, rng_t(j.first(), j.last()-1)^h, k) *= cur_dt / old_dt;
          this->mem->GC[2](i, j, rng_t(k.first(), k.last()-1)^h) *= cur_dt / old_dt;
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
          typename parent_t::bcp_t 
            &bcxl, &bcxr, 
            &bcyl, &bcyr,
            &bczl, &bczr; 
          const rng_t &i, &j, &k; 
        };  

        struct rt_params_t : parent_t::rt_params_t
        {
          typename parent_t::real_t di = 0, dj = 0, dk = 0;
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
	    idx_t<parent_t::n_dims>({args.i, args.j, args.k})
          ),
	  i(args.i), 
          j(args.j), 
          k(args.k),
          stat_field(args.mem->tmp[__FILE__][0][0])
	{
          this->di = p.di;
          this->dj = p.dj;
          this->dk = p.dk;
          this->dijk = {p.di, p.dj, p.dk};
          this->set_bcs(0, args.bcxl, args.bcxr);
	  this->set_bcs(1, args.bcyl, args.bcyr);
	  this->set_bcs(2, args.bczl, args.bczr);
        } 

	public:

	static void alloc(
          typename parent_t::mem_t *mem,
          const int &n_iters
        )   
        {
          // psi
          mem->psi.resize(parent_t::n_eqns);
	  for (int e = 0; e < parent_t::n_eqns; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t(
                parent_t::rng_sclr(mem->grid_size[0]),
                parent_t::rng_sclr(mem->grid_size[1]),
                parent_t::rng_sclr(mem->grid_size[2])
              ))); 

          // Courant field components (Arakawa-C grid)
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t( 
            parent_t::rng_vctr(mem->grid_size[0]),
            parent_t::rng_sclr(mem->grid_size[1]),
            parent_t::rng_sclr(mem->grid_size[2])
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(mem->grid_size[0]),
            parent_t::rng_vctr(mem->grid_size[1]),
            parent_t::rng_sclr(mem->grid_size[2])
          )));
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(mem->grid_size[0]),
            parent_t::rng_sclr(mem->grid_size[1]),
            parent_t::rng_vctr(mem->grid_size[2])
          )));

          // fully third-order accurate mpdata needs also time derivatives of
          // the Courant field
          if (opts::isset(ct_params_t::opts, opts::div_3rd) ||
              opts::isset(ct_params_t::opts, opts::div_3rd_dt))
          {
            // TODO: why for (auto f : {mem->ndt_GC, mem->ndtt_GC}) doesn't work ?
            mem->ndt_GC.push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_vctr(mem->grid_size[0]),
              parent_t::rng_sclr(mem->grid_size[1]),
              parent_t::rng_sclr(mem->grid_size[2])
            )));
            mem->ndt_GC.push_back(mem->old(new typename parent_t::arr_t(
              parent_t::rng_sclr(mem->grid_size[0]),
              parent_t::rng_vctr(mem->grid_size[1]),
              parent_t::rng_sclr(mem->grid_size[2])
            )));
            mem->ndt_GC.push_back(mem->old(new typename parent_t::arr_t(
              parent_t::rng_sclr(mem->grid_size[0]),
              parent_t::rng_sclr(mem->grid_size[1]),
              parent_t::rng_vctr(mem->grid_size[2])
            )));
            
            mem->ndtt_GC.push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_vctr(mem->grid_size[0]),
              parent_t::rng_sclr(mem->grid_size[1]),
              parent_t::rng_sclr(mem->grid_size[2])
            )));
            mem->ndtt_GC.push_back(mem->old(new typename parent_t::arr_t(
              parent_t::rng_sclr(mem->grid_size[0]),
              parent_t::rng_vctr(mem->grid_size[1]),
              parent_t::rng_sclr(mem->grid_size[2])
            )));
            mem->ndtt_GC.push_back(mem->old(new typename parent_t::arr_t(
              parent_t::rng_sclr(mem->grid_size[0]),
              parent_t::rng_sclr(mem->grid_size[1]),
              parent_t::rng_vctr(mem->grid_size[2])
            )));
          }

          // allocate G
          if (opts::isset(ct_params_t::opts, opts::nug))
	    mem->G.reset(mem->old(new typename parent_t::arr_t(
                    parent_t::rng_sclr(mem->grid_size[0]),
                    parent_t::rng_sclr(mem->grid_size[1]),
                    parent_t::rng_sclr(mem->grid_size[2])
            )));

	  // allocate Kahan summation temporary vars
	  if (opts::isset(ct_params_t::opts, opts::khn))
	    for (int n = 0; n < 3; ++n) 
	      mem->khn_tmp.push_back(mem->old(new typename parent_t::arr_t( 
	        parent_t::rng_sclr(mem->grid_size[0]), 
	        parent_t::rng_sclr(mem->grid_size[1]),
	        parent_t::rng_sclr(mem->grid_size[2])
	      )));
          // courant field
          alloc_tmp_sclr(mem, __FILE__, 1);
        }  
        
        // helper method to allocate a temporary space composed of arbitrarily staggered arrays
        static void alloc_tmp_stgr(
          typename parent_t::mem_t *mem,
          const char * __file__,
          const int n_arr,
          const std::vector<std::vector<bool>> &stgr,
          bool srfc = false // allocate only surface data
        )
        {
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          for (int n = 0; n < n_arr; ++n)
          {
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
              stgr[n][0] ? parent_t::rng_vctr(mem->grid_size[0]) : parent_t::rng_sclr(mem->grid_size[0]),
              stgr[n][1] ? parent_t::rng_vctr(mem->grid_size[1]) : parent_t::rng_sclr(mem->grid_size[1]),
              srfc ? rng_t(0, 0) :
                stgr[n][2] ? parent_t::rng_vctr(mem->grid_size[2]) :
                  parent_t::rng_sclr(mem->grid_size[2])
            ))); 
          }
        }
        
        // helper method to allocate a temporary space composed of vector-component arrays
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem,
          const char * __file__
        )
        {
          alloc_tmp_stgr(mem, __file__, 3, {{true, false, false}, {false, true, false}, {false, false, true}});
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem,
          const char * __file__, const int n_arr,
          std::string name = "",
          bool srfc = false // allocate only surface data
        )   
        {   
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());

          if (!name.empty()) mem->avail_tmp[name] = std::make_pair(__file__, mem->tmp[__file__].size() - 1);

          for (int n = 0; n < n_arr; ++n)
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_sclr(mem->grid_size[0]),
              parent_t::rng_sclr(mem->grid_size[1]),
              srfc ? rng_t(0, 0) : parent_t::rng_sclr(mem->grid_size[2])
            )));
        } 
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
