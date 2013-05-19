#include "cloud_common.hpp"

#include <libcloudph++/blk_1m/options.hpp>
#include <libcloudph++/blk_1m/adjustments.hpp>
#include <libcloudph++/blk_1m/forcings_elementwise.hpp>
#include <libcloudph++/blk_1m/forcings_columnwise.hpp>

using namespace advoocat; // TODO: not here?

// @brief a minimalistic kinematic cloud model with bulk microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <
  typename real_t, 
  int n_iters, 
  solvers::inhomo_e inhomo,
  int rhod_th_ix, // dry static energy density divided by c_pd (= dry air density times theta)
  int rhod_rv_ix, // water vapour density
  int rhod_rc_ix, // cloud water density
  int rhod_rr_ix, // rain water density
  int n_eqs = 4
>
class cloud : public cloud_common<real_t, n_iters, inhomo, n_eqs>
{
  using parent_t = cloud_common<real_t, n_iters, inhomo, n_eqs>;

  void condevap()
  {
    auto 
      rhod_th = this->state(rhod_th_ix)(this->ijk),
      rhod_rv = this->state(rhod_rv_ix)(this->ijk),
      rhod_rc = this->state(rhod_rc_ix)(this->ijk),
      rhod_rr = this->state(rhod_rr_ix)(this->ijk);
    auto const
      rhod    = this->rhod(this->ijk);
      
    libcloudphxx::blk_1m::adjustments<real_t>( 
      opts, rhod, rhod_th, rhod_rv, rhod_rc, rhod_rr
    );
  }

  protected:

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    condevap();
    parent_t::hook_ante_loop(nt); // forcings after adjustments
  }

  //
  void update_forcings(typename parent_t::arrvec_t &rhs)
  {
    parent_t::update_forcings(rhs);
 
    auto 
      drhod_rc = rhs.at(rhod_rc_ix),
      drhod_rr = rhs.at(rhod_rr_ix);
    const auto 
      rhod_rc  = this->state(rhod_rc_ix),
      rhod_rr  = this->state(rhod_rr_ix),
      rhod     = this->rhod;

    // element-wise
    {
      const rng_t &i = this->i, &j = this->j;
      libcloudphxx::blk_1m::forcings_elementwise<real_t>(opts, drhod_rc(i,j), drhod_rr(i,j), rhod(i,j), rhod_rc(i,j), rhod_rr(i,j));
    }

    // column-wise
    {
      const rng_t j = this->j;
      for (int i = this->i.first(); i <= this->i.last(); ++i)
	libcloudphxx::blk_1m::forcings_columnwise<real_t>(opts, drhod_rr(i,j), rhod(i,j), rhod_rr(i,j), dz);
    }
  }

  // 
  void hook_post_step()
  {
    condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment
    parent_t::hook_post_step(); // includes the above forcings
    // TODO: shouldn't condevap() be called again here to ensure adjusted field is output?
  }

  real_t dz;
  libcloudphxx::blk_1m::opts<real_t> opts;

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t dz = 0;
    libcloudphxx::blk_1m::opts<real_t> bulk_opts;
  };

  // ctor
  cloud( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    dz(p.dz),
    opts(p.bulk_opts)
  {
    assert(p.dz != 0);
    assert(p.dt != 0);
    opts.dt = p.dt;
  }  
};
