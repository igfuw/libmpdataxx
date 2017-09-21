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
      // 
      template <class ct_params_t > 
      class mpdata_rhs_vip_common : public mpdata_rhs<ct_params_t>
      {
        using ix = typename ct_params_t::ix;

	protected:
	
        using parent_t = mpdata_rhs<ct_params_t>;

	// member fields
        std::array<int, parent_t::n_dims> vip_ixs;
	arrvec_t<typename parent_t::arr_t> &stash, &vip_rhs;
        typename parent_t::real_t eps;

	virtual void fill_stash() = 0;
	virtual void fill_stash_helper(const int d, const int e) final
	{
	  if (ix::vip_den == -1)
	    this->stash[d](this->ijk) = this->state(e)(this->ijk);
	  else if (eps == 0) // this is the default  
          {
            // for those simulations advecting momentum where the division by mass will not cause division by zero
            // (for shallow water simulations it means simulations with no collapsing/inflating shallow water layers)
            this->stash[d](this->ijk) = this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk);
          }
	  else
	  {  
	    this->stash[d](this->ijk) = where(
	      // if
	      this->state(ix::vip_den)(this->ijk) > eps,
	      // then
	      this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk),
	      // else
	      0
	    );
	  }

          assert(std::isfinite(sum(this->stash[d](this->ijk))));
	}

	void extrp(const int d, const int e) // extrapolate velocity field in time to t+1/2
	{                 
	  using namespace arakawa_c;

          const auto beta = this->dt / (2 * this->prev_dt);

          // for dt constant in time we can
          // write the result to stash since we don't need the previous state any more
          // for dt variable in time, however, we have to perform multiple
          // extrapolations per time step and we need to keep the previous state
          // hence the need for another output variable masqueraded as stash[d + offset]
          const auto outd = d + (ct_params_t::var_dt ? ct_params_t::n_dims : 0);

          if (!ct_params_t::var_dt)
          {
	    this->stash[outd](this->ijk) *= -beta;
          }
          else
          {
	    this->stash[outd](this->ijk) = -beta * this->stash[d](this->ijk);
          }

	  if (ix::vip_den == -1) 
	    this->stash[outd](this->ijk) += (1 + beta) * this->state(e)(this->ijk);
	  else if (eps == 0) //this is the default
          {             
            // for those simulations advecting momentum where the division by mass will not cause division by zero
            // (for shallow water simulations it means simulations with no collapsing/inflating shallow water layers)
	    this->stash[outd](this->ijk) += (1 + beta) * (this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk)); 
          }
	  else
	  {
	    this->stash[outd](this->ijk) += where(
	      // if
	      this->state(ix::vip_den)(this->ijk) > eps,
	      // then
	      (1 + beta) * this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk),
	      // else
	      0   
	    );  
	  }

          assert(std::isfinite(sum(this->stash[outd](this->ijk))));
	}   

	arrvec_t<typename parent_t::arr_t>& vips()
        {
          static thread_local arrvec_t<typename parent_t::arr_t> ret;
          for (int d = 0; d < parent_t::n_dims; ++d)
            ret.insert(ret.begin() + d, this->mem->never_delete(&(this->state(vip_ixs[d]))));
          return ret;
        }

	virtual void extrapolate_in_time() = 0;
	virtual void interpolate_in_space() = 0;

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

	void hook_ante_loop(const int nt)
	{
          // fill Courant numbers with zeros so that the divergence test does no harm
          if (this->rank == 0)
            for (int d=0; d < parent_t::n_dims; ++d) this->mem->GC.at(d) = 0; 
          this->mem->barrier();

	  parent_t::hook_ante_loop(nt);
	  
	  // to make extrapolation possible at the first time-step
	  fill_stash();
          vip_rhs_impl_init();
	}

        bool calc_gc() final
        {
	  //extrapolate velocity field in time (t+1/2)
	  extrapolate_in_time();

	  //interpolate from velocity field to courant field (mpdata needs courant numbers from t+1/2)
	  interpolate_in_space();

          // TODO: why???
	  this->mem->barrier();

          return true;
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
	  parent_t::alloc_tmp_sclr(mem, __FILE__, (ct_params_t::var_dt ? 2 : 1) * parent_t::n_dims); // stash
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
