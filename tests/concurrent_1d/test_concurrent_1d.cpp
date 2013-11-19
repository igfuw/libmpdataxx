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

#include <libmpdata++/solvers/adv/mpdata_1d.hpp>
#include <libmpdata++/bcond/cyclic_1d.hpp>

int main()
{
  using namespace libmpdataxx;

  std::cerr << "OpenMP: ";
#if defined(_OPENMP)
  std::cerr << "on" << std::endl;
#else 
  std::cerr << "off" << std::endl;
#endif

  using real_t = long double;
  const int nx = 10, nt = 1000;
   
  // OpenMP
  std::cerr << "OpenMP run" << std::endl;
  {
    using solver_t = solvers::mpdata_1d<real_t>;
    typename solver_t::params_t p;
    p.span = {nx};
    concurr::openmp<solver_t, bcond::cyclic> slv(p);
    slv.advance(nt);
  }

  // Boost.Thread
  std::cerr << "Boost.Thread run" << std::endl;
  {
    using solver_t = solvers::mpdata_1d<real_t>;
    typename solver_t::params_t p;
    p.span = {nx};
    concurr::boost_thread<solver_t, bcond::cyclic> slv(p);
    slv.advance(nt);
  }

  // trheads (i.e. auto)
  std::cerr << "threads run" << std::endl;
  {
    using solver_t = solvers::mpdata_1d<real_t>;
    typename solver_t::params_t p;
    p.span = {nx};
    concurr::threads<solver_t, bcond::cyclic> slv(p);
    slv.advance(nt);
  }

  // serial
  std::cerr << "serial run" << std::endl;
  {
    using solver_t = solvers::mpdata_1d<real_t>;
    typename solver_t::params_t p;
    p.span = {nx};
    concurr::serial<solver_t, bcond::cyclic> slv(p);
    slv.advance(nt);
  }
}
