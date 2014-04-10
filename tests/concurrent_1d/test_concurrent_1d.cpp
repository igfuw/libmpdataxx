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

int main()
{
  using namespace libmpdataxx;

  std::cerr << "OpenMP: ";
#if defined(_OPENMP)
  std::cerr << "on" << std::endl;
#else 
  std::cerr << "off" << std::endl;
#endif

  struct ct_params_t : ct_params_default_t
  { 
    using real_t = long double; 
    enum { n_dims = 1 };
    enum { n_eqs = 1 };
    enum { opts = 0 };
  };

  const int nx = 10, nt = 1000;
   
  // OpenMP
  std::cerr << "OpenMP run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::openmp<solver_t, bcond::cyclic, bcond::cyclic> slv(p);
    slv.advance(nt);
  }

  // Boost.Thread
  std::cerr << "Boost.Thread run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::boost_thread<solver_t, bcond::cyclic, bcond::cyclic> slv(p);
    slv.advance(nt);
  }

  // trheads (i.e. auto)
  std::cerr << "threads run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(p);
    slv.advance(nt);
  }

  // serial
  std::cerr << "serial run" << std::endl;
  {
    using solver_t = solvers::mpdata<ct_params_t>;
    typename solver_t::rt_params_t p;
    p.grid_size = {nx};
    concurr::serial<solver_t, bcond::cyclic, bcond::cyclic> slv(p);
    slv.advance(nt);
  }
}
