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

bool almost_equal(T x, T y, T eps)
{
  double err = (x - y) / (x + y);
  std::cerr << "relative error = " << err << std::endl;
  return err <= eps;
}

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

  decltype(run.advectee()) true_solution0(run.advectee_global().shape());
  decltype(run.advectee()) true_solution1(run.advectee_global().shape());
  true_solution0 = run.advectee_global(0);
  true_solution1 = run.advectee_global(1);

  run.advector(0) = 1.0 * dt/dx;

  run.advance(nt);

  auto L2_error0 = sqrt(sum(pow2(run.advectee_global(0) - true_solution0)) / (np * np)) / time;
  auto L2_error1 = sqrt(sum(pow2(run.advectee_global(1) - true_solution1)) / (np * np)) / time;

  return L2_error0 + L2_error1;
}

int main()
{
#if defined(USE_MPI)
  MPI::Init_thread(MPI_THREAD_MULTIPLE);
#endif

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

  std::map<int, std::string> test_description = {
    {0, "without delayed advection"},
    {1, "with delayed advection of advectee 0"},
    {2, "with delayed advection of advectee 1"},
    {3, "with delayed advection of advectees 1 and 2"}
  };

  for(int i=0;i<4;++i)
    std::cerr << "error value in test " << test_description[i] << " = " << std::setprecision(20) << err[i] << std::endl;

  for(int i=1;i<4;++i){
    std::cerr << "comparing error values in the test " << test_description[i] << " and in the test " << test_description[0] << "..." << std::endl;
    if(!almost_equal(err[i], err[0], 1e-13))
    {
      throw std::runtime_error("delay error");
    }
  }

#if defined(USE_MPI)
  MPI::Finalize();
#endif
}
