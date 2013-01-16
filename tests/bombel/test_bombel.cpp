/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a bombel 
 */

// advection
#include <advoocat/lib.hpp>
#include <advoocat/solver_inhomo.hpp>

// plotting
// (<> should be used instead of "" in normal usage) 
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

// auto-deallocating containers
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

// 
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

enum {tht, prs};
enum {x, y};

template <int n_iters>
class bombel : public inhomo_solver<solvers::mpdata_2d<n_iters, cyclic_2d<x>, cyclic_2d<y>,2>>
{
  void forcings(real_t dt)
  {
  //TODO
  }

  public:

  bombel(int nx, int ny, real_t dt) :
    inhomo_solver<solvers::mpdata_2d<n_iters, cyclic_2d<x>, cyclic_2d<y>, 2>>(nx, ny, dt)
  {
  }
};

int main() 
{
  const int nx = 100, ny = 100, nt = 10, n_out=1;
  const real_t C = .5, dt = 1;
 
  bombel<2> solver(nx, ny, dt);

  // initial condition
  {
    blitz::firstIndex i;
    blitz::firstIndex j;
    solver.state(tht) = real_t(0); 
    solver.state(prs) = real_t(0);
  }
  solver.courant(x) = C;
  solver.courant(y) = C;

  // integration
  for (int t = n_out; t <= nt; t+=n_out)
  {
    solver.solve(n_out);
  }
};
