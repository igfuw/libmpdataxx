/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

//#include <libmpdata++/.../mpdata_3d.hpp>
#include <libmpdata++/solvers/adv/donorcell_3d.hpp>
#include <libmpdata++/bcond/cyclic_3d.hpp> // TODO: needed?
#include <libmpdata++/concurr/threads.hpp>

enum {x, y, z};

int main() 
{
  using namespace libmpdataxx;

  using solver_t = solvers::donorcell_3d<float>;
  typename solver_t::params_t p;

  p.span = {24, 24, 24};

  concurr::threads<
    solver_t,
    bcond::cyclic,
    bcond::cyclic,
    bcond::cyclic
  > slv();
  // TODO: test mpdata
}
