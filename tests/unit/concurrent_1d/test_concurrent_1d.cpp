/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "concurrent_1d/test_concurrent_1d.cpp"
 */

#include <libmpdata++/concurr/openmp.hpp>
#include <libmpdata++/concurr/boost_thread.hpp>
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/concurr/threads.hpp>

#include <libmpdata++/solvers/mpdata.hpp>

#include <libmpdata++/output/gnuplot.hpp>

int main()
{
  using namespace libmpdataxx;

  std::cerr << "OpenMP: ";
#if defined(_OPENMP)
  std::cerr << "on" << std::endl;
#else 
  std::cerr << "off" << std::endl;
#endif

//<listing-1>
  struct ct_params_t : ct_params_default_t
  { 
    using real_t = double; 
    enum { n_dims = 1 }; 
    enum { n_eqns = 1 };
  };
//</listing-1>

  const int nx = 10, nt = 1000;
   
  // OpenMP
  std::cerr << "OpenMP run" << std::endl;
  {
//<listing-2>
    using slv_t = solvers::mpdata<ct_params_t>;
//</listing-2>
//<listing-3>
    using slv_out_t = output::gnuplot<slv_t>;
//</listing-3>
//<listing-4>
    using run_t = concurr::openmp<
      slv_out_t, 
      bcond::cyclic, bcond::cyclic
    >;
//</listing-4>
//<listing-5>
    typename slv_out_t::rt_params_t p;
    p.grid_size = { nx };
    run_t run(p);
//</listing-5>
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }

  // Boost.Thread
  std::cerr << "Boost.Thread run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::boost_thread<solver_t, bcond::cyclic, bcond::cyclic> run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }

  // trheads (i.e. auto)
  std::cerr << "threads run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }

  // serial
  std::cerr << "serial run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::serial<solver_t, bcond::cyclic, bcond::cyclic> run(p);
    run.advectee() = 0;
    run.advector() = 0;
    run.advance(nt);
  }
}
