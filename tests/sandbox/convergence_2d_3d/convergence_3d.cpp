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
#include <libmpdata++/concurr/threads.hpp>

using namespace libmpdataxx;

using T = double;

template <int opts_arg, int opts_iters>
T test(int np)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { n_dims = 3 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };

  typename ct_params_t::real_t
    time = 1.0,
    dx = 1.0 / (np - 1),  
    dy = 1.0 / (np - 1),
    dz = 1.0 / (np - 1),
    dt = 0.1 * dx,
    pi = boost::math::constants::pi<typename ct_params_t::real_t>();

  int nt = time / dt;

  using slv_t = solvers::mpdata<ct_params_t>;

  typename slv_t::rt_params_t p;

  p.n_iters = opts_iters; 
  p.grid_size = {np, np, np};

  concurr::threads<
    slv_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p); 


  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;

  run.advectee() = 2 + sin(2 * pi * i * dx) * sin(2 * pi * j * dy) * sin(2 * pi * k * dz);

  decltype(run.advectee()) true_solution(run.advectee().shape());
  true_solution = run.advectee();

  run.advector(0) = 1.0 * dt/dx;
  run.advector(1) = 2.0 * dt/dy;
  run.advector(2) = 3.0 * dt/dy;

  run.advance(nt);

  auto L2_error = sqrt(sum(pow2(run.advectee() - true_solution)) / (np * np * np)) / time;
  return L2_error;
}

template <int opts_arg, int opts_iters>
int order()
{
    auto err_coarse = test<opts_arg, opts_iters>(17);
    auto err_fine = test<opts_arg, opts_iters>(33);
    return std::round(std::log2(err_coarse / err_fine));
}

int main()
{
  {
    enum { opts = 0 };
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 2) throw std::runtime_error("basic");
  }
  
  {
    enum { opts = opts::iga | opts::fct};
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 2) throw std::runtime_error("default");
  }
  
  {
    enum { opts = opts::tot };
    enum { opts_iters = 3};
    if (order<opts, opts_iters>() != 3) throw std::runtime_error("tot");
  }
  
  {
    enum { opts = opts::abs | opts::tot};
    enum { opts_iters = 3};
    if (order<opts, opts_iters>() != 3) throw std::runtime_error("abs_tot");
  }
  
  {
    enum { opts = opts::iga | opts::tot};
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 3) throw std::runtime_error("iga_tot");
  }
}
