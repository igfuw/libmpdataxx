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

template <int opts_delayed_step_arg>
T test()
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = T;
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
    enum { opts = opts::iga | opts::fct };
    enum { opts_delayed_step = opts_delayed_step_arg };
  };

  const int np = 65;

  typename ct_params_t::real_t
    time = 1.0,
    dx = 1.0 / (np - 1),  
    dt = 0.2 * dx,
    pi = boost::math::constants::pi<typename ct_params_t::real_t>();

  int nt = time / dt;

  using slv_t = solvers::mpdata<ct_params_t>;

  typename slv_t::rt_params_t p;

  p.n_iters = 2; 
  p.grid_size = {np};

  concurr::threads<
    slv_t, 
    bcond::cyclic, bcond::cyclic
  > run(p); 


  blitz::firstIndex i;

  run.advectee(0) = 2 + sin(2 * pi * i * dx);
  run.advectee(1) = 2 + sin(4 * pi * i * dx);

  decltype(run.advectee()) true_solution0(run.advectee().shape());
  decltype(run.advectee()) true_solution1(run.advectee().shape());
  true_solution0 = run.advectee(0);
  true_solution1 = run.advectee(1);

  run.advector(0) = 1.0 * dt/dx;

  run.advance(nt);

  auto L2_error0 = sqrt(sum(pow2(run.advectee(0) - true_solution0)) / (np * np)) / time;
  auto L2_error1 = sqrt(sum(pow2(run.advectee(1) - true_solution1)) / (np * np)) / time;

  return L2_error0 + L2_error1;
}

int main()
{
  double err[4];
  {
    enum { delayed_step = 0 };
    err[0] = test<delayed_step>(); 
  }
  {
    enum { delayed_step = opts::bit(0) };
    err[1] = test<delayed_step>(); 
  }
  {
    enum { delayed_step = opts::bit(1) };
    err[2] = test<delayed_step>(); 
  }
  {
    enum { delayed_step = opts::bit(0) | opts::bit(1) };
    err[3] = test<delayed_step>(); 
  }

  for(int i=1;i<4;++i){
    if(err[i] != err[i-1]) throw std::runtime_error("delay error");
  }
}
