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
#include "convergence_3d_spacetime.hpp"

using namespace libmpdataxx;

template <int opts_arg, int opts_iters>
double test(int np)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 3 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };

  typename ct_params_t::real_t
    pi = boost::math::constants::pi<typename ct_params_t::real_t>(),
    time = 1.0,
    dx = 2.0 * pi / (np - 1),
    dy = 2.0 * pi / (np - 1),  
    dz = 2.0 * pi / (np - 1), 
    dt = 0.5 * dx;

  int nt = time / dt;
  time = nt * dt;

  using slv_t = convergence_3d_spacetime<ct_params_t>;

  typename slv_t::rt_params_t p;
  p.n_iters = opts_iters; 
  p.grid_size = {np, np, np};
  p.di = dx;
  p.dj = dy;
  p.dk = dz;
  p.dt = dt;

  concurr::threads<
    slv_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p); 

  run.advectee() = 8;

  // just to fool the divergence check
  run.advector(0) = 0;
  run.advector(1) = 0;
  run.advector(2) = 0;
  
  blitz::firstIndex i;
  blitz::secondIndex j;
  blitz::thirdIndex k;

  decltype(run.advectee()) true_solution(run.advectee_global().shape());
  decltype(run.advectee()) g_factor_global(true_solution.shape());
  true_solution = (2 + sin(i * dx) * sin(time)) * (2 + sin(j * dy) * sin(time)) * (2 + sin(k * dz) * sin(time));
  
  run.g_factor() = exp(cos(i * dx) + cos(j * dy) + cos(k * dz));
  g_factor_global = exp(cos(i * dx) + cos(j * dy) + cos(k * dz));

  run.advance(nt);

  auto L2_error = sqrt(sum(g_factor_global * pow2(run.advectee_global() - true_solution)))
                  / sqrt(sum(g_factor_global *  pow2(true_solution)));
  return L2_error;
}

template <int opts_arg, int opts_iters>
int order()
{
    auto err_coarse = test<opts_arg, opts_iters>(33);
    auto err_fine = test<opts_arg, opts_iters>(65);
    return std::round(std::log2(err_coarse / err_fine));
}

int main()
{
#if defined(USE_MPI)
  // we will instantiate many solvers, so we have to init mpi manually, 
  // because solvers will not know should they finalize mpi upon destruction
  MPI::Init_thread(MPI_THREAD_SERIALIZED);
#endif
  // mpdata without dfl option set is frist-order accurate for divergent flows
  {
    enum { opts = opts::nug };
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 1) throw std::runtime_error("basic");
  }
  
  // with dfl option we recover second-order accuracy
  {
    enum { opts = opts::nug | opts::dfl };
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 2) throw std::runtime_error("dfl");
  }
  
  // same using the divergence form of antidiffusive velocity
  {
    enum { opts = opts::nug | opts::div_2nd };
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 2) throw std::runtime_error("div_2nd");
  }
  
  // third-order terms for constant coefficients do not result in third-order convergence for this test
  {
    enum { opts = opts::nug | opts::dfl | opts::tot};
    enum { opts_iters = 3};
    if (order<opts, opts_iters>() != 2) throw std::runtime_error("tot");
  }
  
  // using full third-order correction leads to third-order convergence
  {
    enum { opts = opts::nug | opts::div_2nd | opts::div_3rd};
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 3) throw std::runtime_error("div_3rd");
  }

  // same if we use the infinite gauge option
  {
    enum { opts = opts::nug | opts::iga | opts::div_2nd | opts::div_3rd};
    enum { opts_iters = 2};
    if (order<opts, opts_iters>() != 3) throw std::runtime_error("iga_div_3rd");
  }
#if defined(USE_MPI)
  MPI::Finalize();
#endif

}
