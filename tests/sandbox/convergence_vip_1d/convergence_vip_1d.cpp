/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

/* This test checks the convergence order for the 1D Burgers equation u_t + (uu)_x = 0.
 * The initial condition is u_0(x) = 2 + sin(x). The solution is integrated to time pi / 16,
 * well before shock formation, so we should expect maximal convergence order.
 * The primary purpose of this test is to check that the procedures for interpolation/extrapolation
 * and calculation of temporal derivatives (needed in the case of the fully third-order scheme)
 * work correctly in the dynamical solver mpdata_rhs_vip
 */

#include <libmpdata++/concurr/serial.hpp>
#include <boost/math/constants/constants.hpp>
#include <libmpdata++/solvers/mpdata_rhs_vip.hpp>

using namespace libmpdataxx;
using T = double;

const T pi = boost::math::constants::pi<T>();

// The exact solution to our problem u(t, x) satisfies a nonlinear equation u = 2 + sin(x - 2 * u * t).
// The equation is solved to machine precision by a fixed point iteration.
T exact(const T t, const T x)
{
  T u = std::sin(x); // "arbitrary" starting value
  T r = 1;

  while (r > 1e-15)
  {
    const T nu = 2 + std::sin(x - 2 * u * t);
    r = std::abs(nu - u);
    u = nu;
  }

  return u;
}

template<typename ct_params_t>
class solver : public solvers::mpdata_rhs_vip<ct_params_t>
{
  using parent_t = solvers::mpdata_rhs_vip<ct_params_t>;
  using parent_t::parent_t;

  void hook_ante_loop(const typename parent_t::advance_arg_t nt)
  {
    parent_t::hook_ante_loop(nt);

    // filling vip_stash with analytical data for t < t_0
    // seems to be needed to achive the desired convergence rate only for the third-order scheme
    // but done for the second-order scheme as well for consistency
    this->dt_stash.fill(this->dt);
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      this->vip_stash(-1)[0](i) = exact(    -this->dt, i * this->di);
      if (parent_t::div3_mpdata) this->vip_stash(-2)[0](i) = exact(-2 * this->dt, i * this->di);
    }
  }
};

struct err_t
{
  T li, l2;
};

template <int opts_arg, int opts_iters, bool var_dt_arg>
err_t test(int np)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { var_dt = var_dt_arg };
    enum { opts = opts_arg };
    enum { n_dims = 1 };
    enum { n_eqns = 1 };
    enum { rhs_scheme = solvers::trapez };
    struct ix { enum {
      u,
      vip_i=u, vip_den=-1
    }; };
  
    enum { hint_norhs = opts::bit(ix::u) };
  }; 

  using ix = typename ct_params_t::ix;

  using solver_t = solver<ct_params_t>;
  typename solver_t::rt_params_t p;

  p.di = 2.0 * pi / (np - 1);
  p.n_iters = opts_iters;
  p.grid_size = {np};
  
  T time = pi / 16.;
  p.dt = p.di / 16.;
  p.max_courant = 0.25;
  int nt = time / p.dt;
  time = nt * p.dt;

  libmpdataxx::concurr::serial<
    solver_t, 
    bcond::cyclic, bcond::cyclic
  > slv(p);

  {
    blitz::firstIndex i;
    slv.advectee(ix::u) = 2 + sin(i * p.di);
  }

  slv.advance(ct_params_t::var_dt ? time : nt); 

  T li_err = 0.;
  T li_norm = 0.;
  T l2_err = 0;
  T l2_norm = 0.;
  
  for (int i = 0; i < np; ++i)
  {
    const T ex = exact(slv.time(), i * p.di);
    const T nm = slv.advectee(ix::u)(i);

    li_err = std::max(li_err, std::abs(ex - nm));
    li_norm = std::max(li_norm, std::abs(ex));
    
    l2_err += (ex - nm) * (ex - nm);
    l2_norm += ex * ex;
  }

  li_err /= li_norm;

  l2_err /= l2_norm;
  l2_err = std::sqrt(l2_err);

  return {li_err, l2_err};
};

template <int opts_arg, int opts_iters>
bool check(const int expected)
{
  const int n_coarse = 1024;
  const int n_fine = 2 * (n_coarse - 1) + 1;

  // constant time step
  auto err_coarse_cdt = test<opts_arg, opts_iters, false>(n_coarse);
  auto err_fine_cdt = test<opts_arg, opts_iters, false>(n_fine);

  int ord_cdt_l2 = std::round(std::log2(err_coarse_cdt.l2 / err_fine_cdt.l2));
  int ord_cdt_li = std::round(std::log2(err_coarse_cdt.li / err_fine_cdt.li));

  // variable time step
  auto err_coarse_vdt = test<opts_arg, opts_iters, true>(n_coarse);
  auto err_fine_vdt = test<opts_arg, opts_iters, true>(n_fine);

  int ord_vdt_l2 = std::round(std::log2(err_coarse_vdt.l2 / err_fine_vdt.l2));
  int ord_vdt_li = std::round(std::log2(err_coarse_vdt.li / err_fine_vdt.li));

  return ord_cdt_l2 == expected && ord_cdt_li == expected && ord_vdt_l2 == expected && ord_vdt_li == expected;
}

int main()
{
  if (!check<opts::iga | opts::dfl, 2>(2)) throw std::runtime_error("");
  if (!check<opts::iga | opts::div_2nd | opts::div_3rd, 2>(3)) throw std::runtime_error("");
}
