#include <advoocat/solvers/mpdata_2d.hpp>
#include <advoocat/solvers/solver_inhomo.hpp>

#include <libcloudph++/bulk/condevap.hpp>
#include <libcloudph++/bulk/autoconv.hpp>
#include <libcloudph++/bulk/collect.hpp>

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
class cloud : public solvers::inhomo_solver<solvers::mpdata_2d<real_t, n_iters, n_eqs>, inhomo>
{
  using parent_t = solvers::inhomo_solver<solvers::mpdata_2d<real_t, n_iters, n_eqs>, inhomo>;

  void condevap()
  {
    auto 
      rhod_th = this->state(rhod_th_ix)(this->i, this->j),
      rhod_rv = this->state(rhod_rv_ix)(this->i, this->j),
      rhod_rc = this->state(rhod_rc_ix)(this->i, this->j),
      rhod_rr = this->state(rhod_rr_ix)(this->i, this->j);

    libcloudphxx::bulk::condevap<real_t>( // TODO: dt as arg needed?
      rhod, rhod_th, rhod_rv, rhod_rc, rhod_rr
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
      drhod_rc = rhs.at(rhod_rc_ix)(this->i, this->j),
      drhod_rr = rhs.at(rhod_rr_ix)(this->i, this->j);
    const auto 
      rhod_rc = this->state(rhod_rc_ix)(this->i, this->j),
      rhod_rr = this->state(rhod_rr_ix)(this->i, this->j);

    libcloudphxx::bulk::autoconv<real_t>(drhod_rc, drhod_rr, rhod, rhod_rc, rhod_rr);
    libcloudphxx::bulk::collect<real_t>(drhod_rc, drhod_rr, rhod, rhod_rc, rhod_rr);
  }

  // 
  void hook_post_step()
  {
    condevap(); // treat saturation adjustment as post-advection, pre-rhs adjustment
    parent_t::hook_post_step(); // includes the above forcings
    // TODO: shouldn't condevap() be called again here to ensure adjusted field is output?
  }

  // with the size/range of the subdomain
  typename parent_t::arr_t rhod;

  public:

  struct params_t : parent_t::params_t 
  { 
    typename parent_t::arr_t rhod;
  };

  // ctor
  cloud( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p)
  {
    // initialising rhod array with data from the p.rhod profile
    rhod.resize(this->i, this->j);

    assert(
      p.rhod.extent(0) == rhod.extent(1) 
      && "p.rhod profile has wrong dimension (or is not initialised)"
    );

    blitz::firstIndex i;
    blitz::secondIndex j;
    rhod = p.rhod(j);
  }  
};
