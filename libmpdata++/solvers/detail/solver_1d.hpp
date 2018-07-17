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

	protected:

	const rng_t i; //TODO: to be removed

        // generic field used for various statistics (currently Courant number and divergence)
        typename parent_t::arr_t &stat_field; // TODO: should be in solver common but cannot be allocated there ?

        virtual void xchng_sclr(typename parent_t::arr_t &arr, const bool deriv = false) final // for a given array
        {
          this->mem->barrier();
          for (auto &bc : this->bcs[0]) bc->fill_halos_sclr(arr, deriv);
          this->mem->barrier();
        }

        // no pressure solver in 1D but this function needs to be present for dimension independant code,
        // should be only used with cyclic boundary conditions where xchng_pres == xchng_sclr
        virtual void xchng_pres(typename parent_t::arr_t &arr, const idx_t<1>&, const int ext = 0) final
        {
          xchng_sclr(arr);
        }

	void xchng(int e) final
	{
          xchng_sclr(this->mem->psi[e][ this->n[e]]);
	}

        void xchng_vctr_alng(arrvec_t<typename parent_t::arr_t> &arrvec, const bool ad = false, const bool cyclic = false) final
        {
          this->mem->barrier();
          if (!cyclic)
          {
            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_alng(arrvec, ad); 
          }
          else
          {
            for (auto &bc : this->bcs[0]) bc->fill_halos_vctr_alng_cyclic(arrvec, ad);
          }
          this->mem->barrier();
        }

        typename parent_t::real_t courant_number(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          stat_field(this->ijk) = 0.5 * (abs(arrvec[0](i+h) + arrvec[0](i-h)));
          return this->mem->max(this->rank, stat_field(this->ijk));
        }
        
        typename parent_t::real_t max_abs_vctr_div(const arrvec_t<typename parent_t::arr_t> &arrvec) final
        {
          stat_field(this->ijk) = abs((arrvec[0](i+h) - arrvec[0](i-h)));
          return this->mem->max(this->rank, stat_field(this->ijk));
        }

        void scale_gc(const typename parent_t::real_t time,
                      const typename parent_t::real_t cur_dt,
                      const typename parent_t::real_t old_dt) final
        {
          this->mem->GC[0](rng_t(i.first(), i.last()-1)^h) *= cur_dt / old_dt;
          this->xchng_vctr_alng(this->mem->GC);
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
          typename parent_t::real_t di = 0;
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
          i(args.i),
          stat_field(args.mem->tmp[__FILE__][0][0])
	{
          this->di = p.di;
          this->dijk = {p.di};
          this->set_bcs(0, args.bcxl, args.bcxr);
        }

        // memory allocation logic using static methods

        private:

	static void alloc_tmp(
	  typename parent_t::mem_t *mem, 
	  const char * __file__, 
	  const int n_arr,
	  const rng_t rng,
          std::string name = ""
	)
	{
	  mem->tmp[__file__].push_back(new arrvec_t<typename parent_t::arr_t>()); 

          if (!name.empty()) mem->avail_tmp[name] = std::make_pair(__file__, mem->tmp[__file__].size() - 1);

          for (int n = 0; n < n_arr; ++n)
          {
	    mem->tmp[__file__].back().push_back(
              mem->old(new typename parent_t::arr_t( rng ))
            ); 
          }
	}

        public:

	static void alloc(
          typename parent_t::mem_t *mem, 
          const int &n_iters
        ) {
          mem->psi.resize(parent_t::n_eqns);
	  for (int e = 0; e < parent_t::n_eqns; ++e) // equations
	    for (int n = 0; n < n_tlev; ++n) // time levels
	      mem->psi[e].push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_sclr(mem->grid_size[0]))));
    
	  mem->GC.push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_vctr(mem->grid_size[0])))); 

          // fully third-order accurate mpdata needs also time derivatives of
          // the Courant field
          if (opts::isset(ct_params_t::opts, opts::div_3rd) ||
              opts::isset(ct_params_t::opts, opts::div_3rd_dt))
          {
            // TODO: why for (auto f : {mem->ndt_GC, mem->ndtt_GC}) doesn't work ?
	    mem->ndt_GC.push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_vctr(mem->grid_size[0]))));
	    mem->ndtt_GC.push_back(mem->old(new typename parent_t::arr_t(parent_t::rng_vctr(mem->grid_size[0]))));
          }

          if (opts::isset(ct_params_t::opts, opts::nug))
	    mem->G.reset(mem->old(new typename parent_t::arr_t(parent_t::rng_sclr(mem->grid_size[0]))));

          // allocate Kahan summation temporary vars
          if (opts::isset(ct_params_t::opts, opts::khn))
            for (int n = 0; n < 3; ++n) 
              mem->khn_tmp.push_back(mem->old(new typename parent_t::arr_t( 
                parent_t::rng_sclr(mem->grid_size[0])
              )));

          // courant field
          alloc_tmp_sclr(mem, __FILE__, 1);
        } 

        protected:

        // helper method to allocate a vector-component temporary array
        static void alloc_tmp_vctr(
          typename parent_t::mem_t *mem, 
          const char * __file__
        )
        {
          alloc_tmp(mem, __file__, 1, parent_t::rng_vctr(mem->grid_size[0])); // always one-component vectors
        }

        // helper method to allocate n_arr scalar temporary arrays 
        static void alloc_tmp_sclr(
          typename parent_t::mem_t *mem, 
          const char * __file__, const int n_arr
        )
        {
          alloc_tmp(mem, __file__, n_arr, parent_t::rng_sclr(mem->grid_size[0])); 
        }
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
