/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a minimalistic model of a harmonic oscillator
 * (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
 *
 * \image html "../../tests/harmosc/figure_euler_it=1.svg"
 * \image html "../../tests/harmosc/figure_euler_it=1.svg"
 * \image html "../../tests/harmosc/figure_strang_it=2.svg"
 * \image html "../../tests/harmosc/figure_strang_it=2.svg"
 */

#include "coupled_harmosc.hpp"
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx;

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

using real_t = double;

const int nt = 1500;

int main() 
{
  struct ct_params_t
  {
    using real_t = real_t;
    enum { n_dims = 1 };
    enum { n_eqs = 2 };
    enum { opts = 0 };
    enum { rhs_scheme = solvers::rhs_scheme_t::strang };
    struct ix { enum {psi, phi}; };
  };

  using solver_t = output::gnuplot<coupled_harmosc<ct_params_t>>;
  typename solver_t::rt_params_t p; 

  // run-time parameters
  p.span = {1000};
  p.dt = 1;
  p.omega = 2*pi<real_t>() / p.dt / 400;
  p.outfreq = 10;

  using ix = typename ct_params_t::ix;
  p.outvars = {
    {ix::psi, {.name = "psi", .unit = "1"}},
    {ix::phi, {.name = "phi", .unit = "1"}}
  };
  p.gnuplot_command = "plot";

  // instantiation
  concurr::threads<solver_t, bcond::cyclic> run(p);

  // initial condition
  {
    blitz::firstIndex i;
    run.advectee(ix::psi) = pow(
      sin((i+.5) * pi<real_t>() / p.span[0] + pi<real_t>()/3), 
      300
    );
    run.advectee(ix::phi) = real_t(0);
  }
  run.advector() = .5;

  // integration
  run.advance(nt);
}
