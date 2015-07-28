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

    namespace detail
    {
      // 
      template <class ct_params_t > 
      class mpdata_rhs_vip_common : public mpdata_rhs<ct_params_t>
      {
	using parent_t = mpdata_rhs<ct_params_t>;
        using ix = typename ct_params_t::ix;

	protected:

	// member fields
	arrvec_t<typename parent_t::arr_t> &stash;
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
	{                 // (write the result to stash since we don't need previous state any more)
	  using namespace arakawa_c;

	  this->stash[d](this->ijk) /= -2.;

	  if (ix::vip_den == -1) 
	    this->stash[d](this->ijk) += 3./2 * this->state(e)(this->ijk);
	  else if (eps == 0) //this is the default
          {             
            // for those simulations advecting momentum where the division by mass will not cause division by zero
            // (for shallow water simulations it means simulations with no collapsing/inflating shallow water layers)
	    this->stash[d](this->ijk) += 3./2 * (this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk)); 
          }
	  else
	  {
	    this->stash[d](this->ijk) += where(
	      // if
	      this->state(ix::vip_den)(this->ijk) > eps,
	      // then
	      3./2 * this->state(e)(this->ijk) / this->state(ix::vip_den)(this->ijk),
	      // else
	      0   
	    );  
	  }

          assert(std::isfinite(sum(this->stash[d](this->ijk))));
	}   

	virtual void extrapolate_in_time() = 0;
	virtual void interpolate_in_space() = 0;
	virtual void add_absorber() { throw std::logic_error("absorber not yet implemented in 1d & 3d"); }
	virtual void add_relax() { throw std::logic_error("absorber not yet implemented in 1d & 3d"); }
	virtual void normalize_absorber() { throw std::logic_error("absorber not yet implemented in 1d & 3d"); }

	void hook_ante_loop(const int nt)
	{
          // fill Courant numbers with zeros so that the divergence test does no harm
          if (this->rank == 0)
            for (int d=0; d < parent_t::n_dims; ++d) this->mem->GC.at(d) = 0; 
          this->mem->barrier();

	  parent_t::hook_ante_loop(nt);
	  
	  // to make extrapolation possible at the first time-step
	  fill_stash();
	}

	void hook_ante_step()
	{ 
	  //extrapolate velocity field in time (t+1/2)
	  extrapolate_in_time();

	  //interpolate from velocity field to courant field (mpdata needs courant numbers from t+1/2)
	  interpolate_in_space();

          // TODO: why???
	  this->mem->barrier();

	  // filling the stash with data from current velocity field 
	  // (so that in the next time step they can be used for extrapolation in time)
	  fill_stash();

	  // intentionally after stash !!!
	  // (we have to stash data from the current time step before applying any forcings to it)
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) != noab) add_absorber();
	  parent_t::hook_ante_step(); 
	}
	
        void hook_post_step()
	{ 
	  parent_t::hook_post_step(); 
          if (static_cast<vip_vab_t>(ct_params_t::vip_vab) == impl)
          {
            add_relax();
            normalize_absorber();
          }
        }

	public:
	
        struct rt_params_t : parent_t::rt_params_t
        {
          typename parent_t::real_t vip_eps = 0;
        };

	static void alloc(
	  typename parent_t::mem_t *mem, 
	  const rt_params_t &p
	) {
	  // psi[n-1] secret stash for velocity extrapolation in time
	  parent_t::alloc(mem, p);
	  parent_t::alloc_tmp_sclr(mem, p.grid_size, __FILE__, parent_t::n_dims); 
	}
 
        protected:

	// ctor
	mpdata_rhs_vip_common(
	  typename parent_t::ctor_args_t args,
	  const rt_params_t &p
	) : 
	  parent_t(args, p),
	  stash(args.mem->tmp[__FILE__][0]),
          eps(p.vip_eps)
	{}
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
