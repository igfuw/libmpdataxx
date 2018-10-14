/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include "solver.hpp"

using namespace libmpdataxx;

struct err_t
{
  real_t l2_err, li_err;
};

template <int opts_arg, int opts_iters, int diffusive_vel_ord>
err_t test(int np)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = ::real_t;
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
    enum { opts = opts_arg };
    enum { rhs_scheme = solvers::trapez };

    struct ix { enum {
      psi, u,
      vip_i=u, vip_den=-1
    }; };

    enum {hint_norhs = opts::bit(ix::psi) | opts::bit(ix::u)};
  };

  typename ct_params_t::real_t
    pi = boost::math::constants::pi<typename ct_params_t::real_t>(),
    time = 1.0,
    dx = 2.0 * pi / (np - 1),  
    dt = 0.4 * dx,
    mu = 0.001;

  int nt = time / dt;
  time = nt * dt;

  using slv_t = adv_diffusion_solver<ct_params_t, diffusive_vel_ord>;

  typename slv_t::rt_params_t p;
  p.n_iters = opts_iters; 
  p.grid_size = {np};
  p.di = dx;
  p.dt = dt;
  p.mu = mu;

  concurr::serial<
    slv_t, 
    bcond::cyclic, bcond::cyclic
  > run(p); 

  blitz::firstIndex i;

  run.advectee(0) = exact_f{0., mu}(i * dx);
  run.advectee(1) = 0;

  decltype(run.advectee()) true_solution(run.advectee().shape());
  true_solution = exact_f{time, mu}(i * dx);

  run.advance(nt);

  auto L2_error = sqrt(sum(pow2(run.advectee() - true_solution)))
                  / sqrt(sum(pow2(true_solution)));
  
  auto Li_error = max(abs(run.advectee() - true_solution));
  return {L2_error, Li_error};
}

template<int opts, int opts_iters, int diffusive_vel_ord>
void convergence()
{
  std::vector<int> nps = {33, 65, 129, 257, 513, 1025, 2049, 4097};

  std::cout << "Checking convergence for " << opts::opts_string(opts)
            << " MPDATA options with " << opts_iters
            << " MPDATA iterations and " << diffusive_vel_ord
            << " order of diffusive velocities" << std::endl;

  std::cout << "Printing grid size, l2 error, li error, l2 order, li_order" << std::endl;
  err_t err_old = {0.0, 0.0};
  for (auto np : nps)
  {
    auto err = test<opts, opts_iters, diffusive_vel_ord>(np);
    std::cout << np << " " 
              << err.l2_err << " " 
              << err.li_err << " " 
              << std::log2(err_old.l2_err / err.l2_err) << " "
              << std::log2(err_old.li_err / err.li_err) << " "
              << std::endl;
    err_old = err;
  }
}

int main()
{
  {
    enum { opts = opts::iga | opts::dfl};
    enum { opts_iters = 2};
    enum { diffusive_vel_ord = 2};
    convergence<opts, opts_iters, diffusive_vel_ord>();
  }

  {
    enum { opts = opts::iga | opts::dfl | opts::tot};
    enum { opts_iters = 2};
    enum { diffusive_vel_ord = 4};
    convergence<opts, opts_iters, diffusive_vel_ord>();
  }
  
  {
    enum { opts = opts::iga | opts::div_2nd | opts::div_3rd};
    enum { opts_iters = 2};
    enum { diffusive_vel_ord = 4};
    convergence<opts, opts_iters, diffusive_vel_ord>();
  }
}
