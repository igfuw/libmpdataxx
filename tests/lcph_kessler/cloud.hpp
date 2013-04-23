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
  real_t acth;

  protected:

  // deals with initials supersaturation
  void hook_ante_loop(int)
  {
std::cerr << "hook_ante_loop()" << std::endl;
    // TODO
    const rng_t &i = this->i, &j = this->j;
    auto rhod_th   = this->state(rhod_th_ix)(i, j);
    auto rhod_rv   = this->state(rhod_rv_ix)(i, j);
    auto rhod_rc   = this->state(rhod_rc_ix)(i, j);
    auto rhod_rr   = this->state(rhod_rr_ix)(i, j);


    libcloudphxx::bulk::condevap<typename parent_t::arr_t, real_t>(
      rhod_th,
      rhod_rv,
      rhod_rc,
      rhod_rr
    );
  }

  // 
  void hook_post_step()
  {
std::cerr << "hook_post_step()" << std::endl;
    parent_t::hook_post_step();

    auto rhod_th   = this->state(rhod_th_ix);
    auto rhod_rv   = this->state(rhod_rv_ix);
    auto rhod_rc   = this->state(rhod_rc_ix);
    auto rhod_rr   = this->state(rhod_rr_ix);

    rng_t &i = this->i, &j = this->j;
    
//    rhod_rv(i,j) += .1; // TODO: call something from bulk_kessler.hpp
  }

  public:

  struct params_t : parent_t::params_t 
  { 
    real_t acth = .001; // TODO: document it! 
  };

  // ctor
  cloud( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p), 
    acth(p.acth)
  {}  
};
