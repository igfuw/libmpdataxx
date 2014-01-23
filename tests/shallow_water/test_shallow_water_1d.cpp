/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "shallow_water.hpp"
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx; 

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

const int nx = 32, nt = 100;

int main() 
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    enum { n_dims = 1 };
    enum { n_eqs = 2 };
 
    // options
    enum { rhs_scheme = solvers::strang };

    // indices
    struct ix { enum {qx, h, vip_i=qx, vip_den=h}; };

    // hints
    enum { hint_norhs = formulae::opts::bit(ix::h) }; 
  };
  using ix = typename ct_params_t::ix;

  // solver & output choice
  using solver_t = output::gnuplot<
    shallow_water<ct_params_t>
  >;

  // run-time parameters
  solver_t::rt_params_t p; 

  p.span = { nx, };

  p.dt = .05;
  p.di = 1;

  p.outfreq = nt / 25;
  p.outvars = {
    {ix::h, { .name = "h",   .unit = "m" }}, 
  };
  p.gnuplot_with = "lines";

  // instantiation
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> run(p);

  // initial condition
  {
    blitz::firstIndex i;

    run.advectee(ix::h) = 1 - .1 * pow(
      sin((i+.5) * pi<typename ct_params_t::real_t>() / nx), // TODO: assumes dx=1
      32
    );

    run.advectee(ix::qx) = 0;
  }

  // integration
  run.advance(nt); 
};
