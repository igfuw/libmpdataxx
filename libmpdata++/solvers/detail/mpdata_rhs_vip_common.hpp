/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include <libmpdata++/solvers/mpdata_rhs.hpp>

namespace libmpdataxx
{
  namespace solvers
  {
    enum vip_vab_t
    {
      noab,
      expl,
      impl
    };

    const std::map<vip_vab_t, std::string> vab2string {
      {noab, "noab"},
      {expl, "expl"},
      {impl, "impl"},
    };

    struct mpdata_rhs_vip_family_tag {};

    namespace detail
    {
      // override default extrapolation/interpolation in ct_params
      template <class ct_params_t> 
      struct ct_params_vip_default_t : ct_params_t
      {
        // only override if it has the default value, preserve any special options
        enum {sptl_intrp = ct_params_t::sptl_intrp > 0 ? ct_params_t::sptl_intrp : aver2};
        // unconditionally override, temporal extrapolation dosen't have any special options for now
        enum {tmprl_extrp = linear2};
      };

      template <class ct_params_t, int minhalo> 
      class mpdata_rhs_vip_common : public mpdata_rhs<ct_params_vip_default_t<ct_params_t>, minhalo>
      {
        using ix = typename ct_params_t::ix;

	protected:
	
        using parent_t = mpdata_rhs<ct_params_vip_default_t<ct_params_t>, minhalo>;

	// member fields
        std::array<int, parent_t::n_dims> vip_ixs;
	arrvec_t<typename parent_t::arr_t> &stash, &vip_rhs;
        typename parent_t::real_t eps;

        arrvec_t<typename parent_t::arr_t>& vip_stash(const int t_lev)
        {
          // t_lev ==  0 -> output for extrapolation/derivatives
          // t_lev == -1 -> (n-1) state
          // t_lev == -2 -> (n-2) state, only available with div3_mpdata
          
          static thread_local arrvec_t<typename parent_t::arr_t> ret;
          ret.resize(parent_t::n_dims);
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            if (!ct_params_t::var_dt && !parent_t::div3_mpdata)
            {
              assert(t_lev == 0 || t_lev == -1);
              // for dt constant in time we can
              // use the same stash since we don't need the previous state any more
              auto& comp = stash[d];
              ret.replace(ret.begin() + d, this->mem->never_delete(&(comp)));
            }
            else if (!parent_t::div3_mpdata)
            {
              // for dt variable in time, however, we have to perform multiple
              // extrapolations per time step and we need to keep the previous state
              assert(t_lev == 0 || t_lev == -1);
              auto &comp = stash[d - t_lev *  ct_params_t::n_dims];
              ret.replace(ret.begin() + d, this->mem->never_delete(&(comp)));
            }
            else
            {
              // for the fully third-order mpdata we need to keep both the (n-1)
              // and the (n-2) state and juggle them around to avoid array copying
              assert(t_lev == 0 || t_lev == -1 || t_lev == -2);
              auto& comp = (t_lev == 0 ? stash[d] :
                this->timestep % 2 == 0 ? stash[d - t_lev * ct_params_t::n_dims] : 
                  stash[d + (3 + t_lev) * ct_params_t::n_dims]
              );
              ret.replace(ret.begin() + d, this->mem->never_delete(&(comp)));
            }
          }
          return ret;
        }

        virtual void fill_stash_helper(const int d) final
	{
          // for third-order mpdata saving to t_lev == -2 so that it becomes -1 at the next time step
          int save_t_lev = parent_t::div3_mpdata ? -2 : -1;
	  if (ix::vip_den == -1)
	    vip_stash(save_t_lev)[d](this->ijk) = vips()[d](this->ijk);
	  else if (eps == 0) // this is the default  
          {
            // for those simulations advecting momentum where the division by mass will not cause division by zero
            // (for shallow water simulations it means simulations with no collapsing/inflating shallow water layers)
            vip_stash(save_t_lev)[d](this->ijk) = vips()[d](this->ijk) / this->state(ix::vip_den)(this->ijk);
          }
	  else
	  {  
	    vip_stash(save_t_lev)[d](this->ijk) = where(
	      // if
	      this->state(ix::vip_den)(this->ijk) > eps,
	      // then
	      vips()[d](this->ijk) / this->state(ix::vip_den)(this->ijk),
	      // else
	      0
	    );
	  }

          assert(std::isfinite(sum(vip_stash(save_t_lev)[d](this->ijk))));
	}

        void fill_stash() 
        {
          for (int d = 0; d < ct_params_t::n_dims; ++d) fill_stash_helper(d);
        }

	void extrp(const int d, const int e) // extrapolate velocity field in time to t+1/2
	{                 
	  using namespace arakawa_c;

          const auto beta = this->dt_stash[0] > 0 ? this->dt / (2 * this->dt_stash[0]) : 0;

          if (!ct_params_t::var_dt && !parent_t::div3_mpdata)
          {
	    this->vip_stash(0)[d](this->ijk) *= -beta;
          }
          else
          {
	    this->vip_stash(0)[d](this->ijk) = -beta * this->vip_stash(-1)[d](this->ijk);
          }

	  if (ix::vip_den == -1) 
	    this->vip_stash(0)[d](this->ijk) += (1 + beta) * this->state(e)(this->ijk);
	  else if (eps == 0) //this is the default
          {             
            // for those simulations advecting momentum where the division by mass will not cause division by zero
            // (for shallow water simulations it means simulations with no collapsing/inflating shallow water layers)
	    this->vip_stash(0)[d](this->ijk) += (1 + beta) * (this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk)); 
          }
	  else
	  {
	    this->vip_stash(0)[d](this->ijk) += where(
	      // if
	      this->state(ix::vip_den)(this->ijk) > eps,
	      // then
	      (1 + beta) * this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk),
	      // else
	      0   
	    );  
	  }

          assert(std::isfinite(sum(this->vip_stash(0)[d](this->ijk))));
	}   

	arrvec_t<typename parent_t::arr_t>& vips()
        {
          static thread_local arrvec_t<typename parent_t::arr_t> ret;
          ret.resize(parent_t::n_dims);
          for (int d = 0; d < parent_t::n_dims; ++d)
            ret.replace(ret.begin() + d, this->mem->never_delete(&(this->state(vip_ixs[d]))));
          return ret;
        }

	virtual void extrapolate_in_time() = 0;
	virtual void interpolate_in_space(arrvec_t<typename parent_t::arr_t> &dst,
                                          const arrvec_t<typename parent_t::arr_t> &src) = 0;

        virtual void vip_rhs_impl_init()
        {
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl)
            {
              vip_rhs[d](this->ijk) = - 
                (*this->mem->vab_coeff)(this->ijk) * (vips()[d](this->ijk) - this->mem->vab_relax[d](this->ijk));
            }
            else
            {
              vip_rhs[d](this->ijk) = 0.0;
            }
          }
        }
        
        virtual void vip_rhs_expl_calc()
        {
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == expl)
          {
            for (int d = 0; d < parent_t::n_dims; ++d)
            {
              // factor of 2 because it is multiplied by 0.5 * dt in vip_rhs_apply
              vip_rhs[d](this->ijk) += -2 * 
                (*this->mem->vab_coeff)(this->ijk) * (vips()[d](this->ijk) - this->mem->vab_relax[d](this->ijk));
            }
          }
        }
        
        virtual void vip_rhs_impl_fnlz()
        {
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl)
          {
            for (int d = 0; d < parent_t::n_dims; ++d)
            {
              vip_rhs[d](this->ijk) = -vips()[d](this->ijk);
            }
            add_relax();
            normalize_vip(vips());
            for (int d = 0; d < parent_t::n_dims; ++d)
            {
              vip_rhs[d](this->ijk) += vips()[d](this->ijk);
              vip_rhs[d](this->ijk) /= (0.5 * this->dt);
            }
          }
        }
        
        void vip_rhs_apply()
        {    
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            vips()[d](this->ijk) += 0.5 * this->dt * vip_rhs[d](this->ijk);
            vip_rhs[d](this->ijk) = 0;
          }
        }

        virtual void normalize_vip(const arrvec_t<typename parent_t::arr_t> &v)
        {
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl)
          {
            for (int d = 0; d < parent_t::n_dims; ++d)
            {
              v[d](this->ijk) /= (1.0 + 0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk));
            }
          }
        }
        
        void add_relax()
        {
          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            this->vips()[d](this->ijk) +=
              0.5 * this->dt * (*this->mem->vab_coeff)(this->ijk) * this->mem->vab_relax[d](this->ijk);
          }
        }

	void hook_ante_loop(const typename parent_t::advance_arg_t nt)
	{
          // fill Courant numbers with zeros so that the divergence test does no harm
          if (this->rank == 0)
            for (int d=0; d < parent_t::n_dims; ++d) this->mem->GC.at(d) = 0; 
          this->mem->barrier();

          for (int d = 0; d < parent_t::n_dims; ++d)
          {
            // these arrays should be multiplied by zeros if used before
            // they are available but to avoid problems with NaNs we fill them with zeros anyway
            vip_stash(-1)[d](this->ijk) = 0;
            if (parent_t::div3_mpdata) vip_stash(-2)[d](this->ijk) = 0;
          }

	  parent_t::hook_ante_loop(nt);
	  
          vip_rhs_impl_init();
	}

        bool calc_gc()
        {
	  //extrapolate velocity field in time (t+1/2)
	  extrapolate_in_time();

	  //interpolate from velocity field to courant field (mpdata needs courant numbers from t+1/2)
	  interpolate_in_space(this->mem->GC, vip_stash(0));

          // TODO: why???
	  this->mem->barrier();

          return true;
        }
        
        void calc_ndt_gc() final
        {
          if (parent_t::div3_mpdata)
          {
            auto ex = this->halo - 1;
            if (this->dt_stash[0] > 0)
            {
              for (int d = 0; d < parent_t::n_dims; ++d)
              {
                vip_stash(0)[d](this->ijk) = this->dt * (vips()[d](this->ijk) - vip_stash(-1)[d](this->ijk)) / this->dt_stash[0];
                this->xchng_pres(this->vip_stash(0)[d], this->ijk, ex);
              }
              interpolate_in_space(this->mem->ndt_GC, vip_stash(0));
              this->mem->barrier();
            }
            
            if (this->dt_stash[0] > 0 && this->dt_stash[1] > 0)
            {
              for (int d = 0; d < parent_t::n_dims; ++d)
              {
                vip_stash(0)[d](this->ijk) = this->dt * this->dt * (
                                                  this->dt_stash[1] * vips()[d](this->ijk)
                                                - (this->dt_stash[1] + this->dt_stash[0]) * vip_stash(-1)[d](this->ijk)
                                                + this->dt_stash[0] * vip_stash(-2)[d](this->ijk)
                                             ) / ( this->dt_stash[0] * this->dt_stash[1]
                                                 * 0.5 * (this->dt_stash[0] + this->dt_stash[1])
                                             );
                this->xchng_pres(this->vip_stash(0)[d], this->ijk, ex);
              }
              interpolate_in_space(this->mem->ndtt_GC, vip_stash(0));
              this->mem->barrier();
            }
          }
        }

	void hook_ante_step()
	{ 
	  // filling the stash with data from current velocity field 
	  // (so that in the next time step they can be used for extrapolation in time)
	  fill_stash();

	  // intentionally after stash !!!
	  // (we have to stash data from the current time step before applying any forcings to it)
          vip_rhs_expl_calc();
          // finish calculating velocity forces before moving on
          this->mem->barrier();

	  parent_t::hook_ante_step(); 
          vip_rhs_apply();
	}
	
        void hook_post_step()
	{ 
	  parent_t::hook_post_step(); 
          vip_rhs_impl_fnlz();
        }

	public:
	
        struct rt_params_t : parent_t::rt_params_t
        {
          typename parent_t::real_t vip_eps = 0;
        };

	static void alloc(
	  typename parent_t::mem_t *mem, 
          const int &n_iters
	) {
	  parent_t::alloc(mem, n_iters);
	  parent_t::alloc_tmp_sclr(mem, __FILE__,
            (parent_t::div3_mpdata ? 3 : ct_params_t::var_dt ? 2 : 1) * parent_t::n_dims); // stash
	  parent_t::alloc_tmp_sclr(mem, __FILE__, parent_t::n_dims); // vip_rhs
	}
 
        protected:

	// ctor
	mpdata_rhs_vip_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) : 
	  parent_t(args, p),
	  stash(args.mem->tmp[__FILE__][0]),
	  vip_rhs(args.mem->tmp[__FILE__][1]),
          eps(p.vip_eps)
	{}
      };
    } // namespace detail
  } // namespace solvers
} // namespace libmpdataxx
