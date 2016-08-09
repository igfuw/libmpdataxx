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

	virtual void xchng_sclr(typename parent_t::arr_t &arr,
	               const idx_t<3> &range_ijk,
                       const int ext = 0,
                       const bool deriv = false
        ) final // for a given array
	{
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_sclr(arr, range_ijk[1]^ext, range_ijk[2]^ext, deriv);
	  for (auto &bc : this->bcs[1]) bc->fill_halos_sclr(arr, range_ijk[2]^ext, range_ijk[0]^ext, deriv);
	  for (auto &bc : this->bcs[2]) bc->fill_halos_sclr(arr, range_ijk[0]^ext, range_ijk[1]^ext, deriv);
          this->mem->barrier();
	}
	void xchng(int e) final
	{
	  this->xchng_sclr(this->mem->psi[e][ this->n[e]], this->ijk, this->halo);
	}

        virtual void xchng_vctr_alng(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_alng(arrvec, j, k); 
          for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_alng(arrvec, k, i); 
          for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_alng(arrvec, i, j);
          this->mem->barrier();
        }

        virtual void xchng_vctr_nrml(
          const arrvec_t<typename parent_t::arr_t> &arrvec,
	  const idx_t<3> &range_ijk,
          const int ext = 0
        ) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_nrml(arrvec[0], range_ijk[2]^ext^1, range_ijk[0]^ext^h);
          for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_nrml(arrvec[0], range_ijk[0]^ext^h, range_ijk[1]^ext^1);

          for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_nrml(arrvec[1], range_ijk[1]^ext^h, range_ijk[2]^ext^1);
          for (auto &bc : this->bcs[2]) bc->fill_halos_vctr_nrml(arrvec[1], range_ijk[0]^ext^1, range_ijk[1]^ext^h);
   
          for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_nrml(arrvec[2], range_ijk[1]^ext^1, range_ijk[2]^ext^h);
          for (auto &bc : this->bcs[1]) bc->fill_halos_vctr_nrml(arrvec[2], range_ijk[2]^ext^h, range_ijk[0]^ext^1);
          this->mem->barrier();
        }

        virtual void xchng_pres(
	  const typename parent_t::arr_t &arr,
	  const idx_t<3> &range_ijk
        ) final
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_pres(arr, range_ijk[1], range_ijk[2]);
          for (auto &bc : this->bcs[1]) bc->fill_halos_pres(arr, range_ijk[2], range_ijk[0]);
          for (auto &bc : this->bcs[2]) bc->fill_halos_pres(arr, range_ijk[0], range_ijk[1]);
          this->mem->barrier();
        }

        virtual void set_edges(
          const arrvec_t<typename parent_t::arr_t> &av,
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
        
        void hook_ante_loop(const int nt) // TODO: this nt conflicts in fact with multiple-advance()-call logic!
        {
          parent_t::hook_ante_loop(nt);
	  
          xchng_vctr_alng(this->mem->GC);
 
          // sanity check for non-divergence of the initial Courant number field
          // TODO: same in 1D
          if (!opts::isset(ct_params_t::opts, opts::dfl))
          {
            typename ct_params_t::real_t max_abs_div = max(abs(
              (
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
              ) / formulae::G<ct_params_t::opts, 0>(*this->mem->G, i, j, k)
            ));
	    if (max_abs_div > this->max_abs_div_eps) 
	      throw std::runtime_error("initial advector field is divergent");
          }
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
          k(args.k)
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
        }  
        
        // helper method to allocate a temporary space composed of vector-component arrays
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem,
          const char * __file__
        )
        {
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_vctr(mem->grid_size[0]),
            parent_t::rng_sclr(mem->grid_size[1]),
            parent_t::rng_sclr(mem->grid_size[2])
          ))); 
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
            parent_t::rng_sclr(mem->grid_size[0]),
            parent_t::rng_vctr(mem->grid_size[1]),
            parent_t::rng_sclr(mem->grid_size[2])
          ))); 
          mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t(
             parent_t::rng_sclr(mem->grid_size[0]),
             parent_t::rng_sclr(mem->grid_size[1]),
             parent_t::rng_vctr(mem->grid_size[2])
          ))); 
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem,
          const char * __file__, const int n_arr,
          std::string name = ""
        )   
        {   
          mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>());

          if (!name.empty()) mem->avail_tmp[name] = std::make_pair(__file__, mem->tmp[__file__].size() - 1);

          for (int n = 0; n < n_arr; ++n)
            mem->tmp[__file__].back().push_back(mem->old(new typename parent_t::arr_t( 
              parent_t::rng_sclr(mem->grid_size[0]),
              parent_t::rng_sclr(mem->grid_size[1]),
              parent_t::rng_sclr(mem->grid_size[2])
            )));
        } 
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
