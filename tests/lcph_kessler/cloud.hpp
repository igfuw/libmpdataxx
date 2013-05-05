#include <advoocat/solvers/mpdata_2d.hpp>
#include <advoocat/solvers/solver_inhomo.hpp>

#include <libcloudph++/bulk/condevap.hpp>
#include <libcloudph++/bulk/autoconv.hpp>
#include <libcloudph++/bulk/collectn.hpp>
#include <libcloudph++/bulk/sediment.hpp>

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
      rhod_th = this->state(rhod_th_ix)(this->ijk),
      rhod_rv = this->state(rhod_rv_ix)(this->ijk),
      rhod_rc = this->state(rhod_rc_ix)(this->ijk),
      rhod_rr = this->state(rhod_rr_ix)(this->ijk);
    auto const
      rhod    = this->rhod(this->ijk);
      
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
      drhod_rc = rhs.at(rhod_rc_ix),
      drhod_rr = rhs.at(rhod_rr_ix);
    const auto 
      rhod_rc  = this->state(rhod_rc_ix),
      rhod_rr  = this->state(rhod_rr_ix),
      rhod     = this->rhod;

    // element-wise
    {
      const rng_t &i = this->i, &j = this->j;
      libcloudphxx::bulk::autoconv<real_t>(drhod_rc(i,j), drhod_rr(i,j), rhod(i,j), rhod_rc(i,j), rhod_rr(i,j));
      libcloudphxx::bulk::collectn<real_t>(drhod_rc(i,j), drhod_rr(i,j), rhod(i,j), rhod_rc(i,j), rhod_rr(i,j));
    }

    // column-wise
    {
      const rng_t j = this->j;
      for (int i = this->i.first(); i <= this->i.last(); ++i)
	libcloudphxx::bulk::sediment<real_t>(drhod_rr(i,j), rhod(i,j), rhod_rr(i,j), dz);
    }
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
  real_t dz;

  public:

  struct params_t : parent_t::params_t 
  { 
    std::vector<real_t> rhod; // profile
    real_t dz = 0;
  };

  // ctor
  cloud( 
    typename parent_t::ctor_args_t args, 
    const params_t &p
  ) : 
    parent_t(args, p),
    dz(p.dz),
    rhod(args.mem->tmp[__FILE__][0][0])
  {
    assert(dz != 0);
    assert(p.rhod.size() == this->j.last()+1);
    // initialising rhod array columnwise with data from the p.rhod profile
    for (int i = this->i.first(); i <= this->i.last(); ++i)
      for (int j = this->j.first(); j <= this->j.last(); ++j)
	rhod(i, j) = p.rhod[j];
  }  

  static void alloc(typename parent_t::mem_t *mem, const int nx, const int ny)
  {
    using namespace advoocat::arakawa_c;
    parent_t::alloc(mem, nx, ny);
    const rng_t i(0, nx-1), j(0, ny-1);
    mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>());
    mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t(i^parent_t::halo, j^parent_t::halo)); // rhod
  }
};
