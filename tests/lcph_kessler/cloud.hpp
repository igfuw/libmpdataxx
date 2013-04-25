#include <advoocat/solvers/mpdata_2d.hpp>
#include <libcloudph++/bulk/condevap.hpp>

using namespace advoocat; // TODO: not here?

// @brief a minimalistic kinematic cloud model with bulk microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <typename real_t, int n_iters, 
  int rhod_th_ix, // dry static energy density divided by c_pd (= dry air density times theta)
  int rhod_rv_ix, // water vapour density
  int rhod_rc_ix, // cloud water density
  int rhod_rr_ix, // rain water density
  int n_eqs = 4
>
class cloud : public solvers::mpdata_2d<real_t, n_iters, n_eqs>
{
  using parent_t = solvers::mpdata_2d<real_t, n_iters, n_eqs>;

  void condevap()
  {
    const rng_t 
      &i = this->i, 
      &j = this->j;

    auto 
      rhod_th   = this->state(rhod_th_ix)(i, j),
      rhod_rv   = this->state(rhod_rv_ix)(i, j),
      rhod_rc   = this->state(rhod_rc_ix)(i, j),
      rhod_rr   = this->state(rhod_rr_ix)(i, j);

    libcloudphxx::bulk::condevap<typename parent_t::arr_t, real_t>(
      rhod, rhod_th, rhod_rv, rhod_rc, rhod_rr
    );
  }
 
  protected:

  // deals with initials supersaturation
  void hook_ante_loop(int nt)
  {
    parent_t::hook_ante_loop(nt); // having it here causes the non-adjusted state to be output @ t=0
    condevap();
  }

  // 
  void hook_post_step()
  {
    parent_t::hook_post_step();
    condevap();
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
