/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include <libmpdata++/solvers/shallow_water.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx; 

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

const int nx = 32, ny = 32, nt = 100;

int main() 
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqns = 3 };

    // options
    enum { rhs_scheme = solvers::trapez };
    
    // indices
    struct ix { enum {qx, qy, h, vip_i=qx, vip_j=qy, vip_den=h}; };

    // hints
    enum { hint_norhs = opts::bit(ix::h) }; 
  };
  using ix = typename ct_params_t::ix;

  // solver & output choice
  using solver_t = output::gnuplot<
    solvers::shallow_water<ct_params_t>
  >;

  // run-time parameters
  solver_t::rt_params_t p; 

  p.grid_size = { nx, ny };

  p.dt = .05;
  p.di = 1;
  p.dj = 1;

  p.outfreq = nt / 25;
  p.outvars = {
    {ix::h, { "h",   "m" }}, 
  };
  p.gnuplot_with = "lines";
  p.gnuplot_output = "figure_%s_%04d.svg";
  p.gnuplot_zrange = p.gnuplot_cbrange = "[.85:1.1]";

  // instantiation
  concurr::threads<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p);

  // initial condition
  {
    blitz::firstIndex i;
    blitz::secondIndex j;

    run.advectee(ix::h) = 1 - .1 * pow(
      sin(i * pi<typename ct_params_t::real_t>() / (nx-1)) * // TODO: assumes dx=dy=1
      sin(j * pi<typename ct_params_t::real_t>() / (ny-1)), 
      32
    );

    run.advectee(ix::qx) = 0;
    run.advectee(ix::qy) = 0;
  }

  // integration
  run.advance(nt); 
};
