/** 
 * @file
 * @example concurrent_1d/test_concurrent_1d.cpp
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "concurrent_1d/test_concurrent_1d.cpp"
 */

#include "advoocat/openmp.hpp"
#include "advoocat/solvers/mpdata_1d.hpp"
#include "advoocat/bcond/cyclic_1d.hpp"
#include "advoocat/sharedmem.hpp"

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
  using mem_t = sharedmem_1d<1, real_t>;
  int nx = 10;
  const int n_iters = 1;
   
  // OpenMP
  {
    openmp<solvers::mpdata_1d<n_iters, bcond::cyclic_1d<real_t>, mem_t>> slv(nx);
    slv.advance(10);
  }

  // Boost.Thread
  {
  }

  // C++11 Thread
  {
  }
}
