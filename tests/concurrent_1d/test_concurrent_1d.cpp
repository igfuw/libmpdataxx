/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "concurrent_1d/test_concurrent_1d.cpp"
 */

#include "advoocat/concurr/openmp.hpp"
#include "advoocat/concurr/boost_thread.hpp"
#include "advoocat/solvers/mpdata_1d.hpp"
#include "advoocat/bcond/cyclic_1d.hpp"

int main()
{
  using namespace advoocat;

  std::cerr << "OpenMP: ";
#if defined(_OPENMP)
  std::cerr << "on" << std::endl;
#else 
  std::cerr << "off" << std::endl;
#endif

  using real_t = long double;
  int nx = 10;
  const int n_iters = 3;
   
  // OpenMP
  std::cerr << "OpenMP run" << std::endl;
  {
    concurr::openmp<
      solvers::mpdata_1d<real_t, n_iters>,
      bcond::cyclic
    > slv(nx);
    slv.advance(1000);
  }

  // Boost.Thread
  std::cerr << "Boost.Thread run" << std::endl;
  {
    concurr::boost_thread<
      solvers::mpdata_1d<real_t, n_iters>,
      bcond::cyclic
    > slv(nx);
    slv.advance(1000);
  }
}
