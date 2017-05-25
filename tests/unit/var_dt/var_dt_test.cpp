/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include "var_dt_test.hpp"

using namespace libmpdataxx;
using T = double;

template <int opts_arg, int opts_iters, bool var_dt_arg>
T test(int np, T cfl)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { var_dt = var_dt_arg };
    enum { n_dims = 1 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };

  T pi = boost::math::constants::pi<typename ct_params_t::real_t>(),
    tau = 256,
    u_0 = 1,
    time = tau * sqrt(2.) / 2,
    dx = 2.0 * pi / (np - 1),
    dt = cfl * dx;
  int nt = time / dt;

  using slv_t = var_dt_test<ct_params_t>;
  typename slv_t::rt_params_t p;

  p.n_iters = opts_iters; 
  p.grid_size = {np};
  p.di = dx;
  p.dt = dt;
  p.max_courant = cfl;

  p.u_0 = u_0;
  p.tau = tau;
  p.pi = pi;

  concurr::serial<
    slv_t, 
    bcond::cyclic, bcond::cyclic
  > run(p); 

  run.advectee() = 2. + sin(blitz::tensor::i * dx);

  // hack to make initial cfl different than chosen so that update_gc triggers at least once
  // TODO: not
  //run.advector(0) = cfl * 0.5;
  run.advance(var_dt_arg ? time : nt);

  decltype(run.advectee()) true_solution(run.advectee().shape());
  true_solution.reindexSelf(run.advectee().lbound());
  true_solution = 2 + sin(blitz::tensor::i * dx - sin(pi * u_0 * run.time() / tau) * tau / pi);

  auto L2 = sqrt(sum(pow2(true_solution - run.advectee())) / sum(pow2(true_solution)));
  return L2;
}

template <int opts_arg, int opts_iters>
bool check(T cfl)
{
    // constant time step
    auto err_coarse_cdt = test<opts_arg, opts_iters, false>(33, cfl);
    auto err_fine_cdt = test<opts_arg, opts_iters, false>(65, cfl);
    int ord_cdt = std::round(std::log2(err_coarse_cdt / err_fine_cdt));

    // variable time step
    auto err_coarse_vdt = test<opts_arg, opts_iters, true>(33, cfl);
    auto err_fine_vdt = test<opts_arg, opts_iters, true>(65, cfl);
    int ord_vdt = std::round(std::log2(err_coarse_vdt / err_fine_vdt));

    return ord_cdt == 2 && ord_vdt == 2 && err_fine_vdt < err_fine_cdt;
}

int main()
{
    enum { opts = 0 };
    enum { opts_iters = 2};

    if (!check<opts, opts_iters>(0.1)) throw std::runtime_error("0.1");
    if (!check<opts, opts_iters>(0.5)) throw std::runtime_error("0.5");
    if (!check<opts, opts_iters>(0.9)) throw std::runtime_error("0.9");
}
