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
    namespace detail
    {
      // 
      template <class ct_params_t > 
      class mpdata_rhs_vip_common : public mpdata_rhs<ct_params_t>
      {
	using parent_t = mpdata_rhs<ct_params_t>;

	protected:

	// member fields
	arrvec_t<typename parent_t::arr_t> &stash;
        bool initial_h_non_zero = false;

	// ctor
	mpdata_rhs_vip_common(
	  typename parent_t::ctor_args_t args,
	  const typename parent_t::rt_params_t &p
	) : 
	  parent_t(args, p),
	  stash(args.mem->tmp[__FILE__][0])
	{}

	virtual void fill_stash() = 0;
	virtual void extrapolate_in_time() = 0;
	virtual void interpolate_in_space() = 0;

	void hook_ante_loop(const int nt)
	{
	  parent_t::hook_ante_loop(nt);
	  
          // set-up initial_h_non_zero
          initial_h_non_zero = min(this->psi_n(ct_params_t::ix::vip_den)(this->ijk)) > 0;

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
	  parent_t::hook_ante_step(); 
	}

	public:
	
	static void alloc(
	  typename parent_t::mem_t *mem, 
	  const typename parent_t::rt_params_t &p
	) {
	  // psi[n-1] secret stash for velocity extrapolation in time
	  parent_t::alloc(mem, p);
	  parent_t::alloc_tmp_sclr(mem, p.span, __FILE__, parent_t::n_dims); 
	}
      };
    }; // namespace detail
  }; // namespace solvers
}; // namespace libmpdataxx
