/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

/** 
 Solver for the constant coefficients advection-diffusion equation
 using the diffusion as advection formulation and utilising
 the prognosed velocity libmpdata++ solver
 */

#pragma once
#include <libmpdata++/solvers/mpdata_rhs_vip.hpp>

using real_t = double;

struct exact_f
{
  real_t t, mu;

  real_t operator()(const real_t x) const
  {
    return 2 + std::sin(x - t) * std::exp(-mu * t);

  }
  BZ_DECLARE_FUNCTOR(exact_f);
};

template <int ord, class arg_t>
inline auto diffusive_vel(
  const arg_t &psi,
  const libmpdataxx::rng_t &i,
  const real_t dx,
  const real_t mu,
  typename std::enable_if<ord == 2>::type* = 0
)
{
  return blitz::safeToReturn(
    -mu * (psi(i + 1) - psi(i - 1)) / (2 * dx * psi(i))
  );
}

template <int ord, class arg_t>
inline auto diffusive_vel(
  const arg_t &psi,
  const libmpdataxx::rng_t &i,
  const real_t dx,
  const real_t mu,
  typename std::enable_if<ord == 4>::type* = 0
)
{
  return blitz::safeToReturn(
    -mu * (-psi(i + 2) + 8 * psi(i + 1) - 8 * psi(i - 1) + psi(i - 2)) / (12 * dx * psi(i))
  );
}


template <class ct_params_t, int diffusive_vel_ord>
class adv_diffusion_solver : public libmpdataxx::solvers::mpdata_rhs_vip<ct_params_t>
{
  using parent_t = libmpdataxx::solvers::mpdata_rhs_vip<ct_params_t>;

  static_assert(diffusive_vel_ord == 2 || diffusive_vel_ord == 4, "unimplemented order");
  
  public:

  using real_t = typename ct_params_t::real_t;
  using ix = typename ct_params_t::ix;

  protected:
  
  real_t mu;

  void hook_ante_loop(const typename parent_t::advance_arg_t nt)
  {
    parent_t::hook_ante_loop(nt);
   
    // convenience
    const auto &i = this->i;
    auto &tmp = this->state(ix::u);
 
    // filling inital stash data using the analytical solution
    // at times t=-dt and t=-2*dt
    
    this->dt_stash.fill(this->dt);

    for (int ii = i.first(); ii <= i.last(); ++ii)
    {
      tmp(ii) = exact_f{-this->dt, mu}(ii * this->di);
    }
    this->xchng_sclr(tmp);
    this->vip_stash(-1)[0](i) = 1 + diffusive_vel<diffusive_vel_ord>(tmp, i, this->di, mu);

    if (parent_t::div3_mpdata)
    {

      for (int ii = i.first(); ii <= i.last(); ++ii)
      {
        tmp(ii) = exact_f{-2 * this->dt, mu}(ii * this->di);
      }
      this->xchng_sclr(tmp);
      this->vip_stash(-2)[0](i) = 1 + diffusive_vel<diffusive_vel_ord>(tmp, i, this->di, mu);
    }
  }
  
  bool calc_gc()
  {
    auto &psi = this->state(ix::psi);
    const auto &i = this->i;

    this->xchng_sclr(psi);
    this->vips()[0](i) = 1 + diffusive_vel<diffusive_vel_ord>(psi, i, this->di, mu);

    return parent_t::calc_gc();
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    real_t mu;
  };

  adv_diffusion_solver(
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) :
  parent_t(args, p),
  mu(p.mu)
  {}
};
