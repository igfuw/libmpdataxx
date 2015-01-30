/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a minimalistic model of a harmonic oscillator
 * (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
 */

#include "coupled_harmosc.hpp"
#include "coupled_harmosc_stats.hpp"
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
#include <boost/math/constants/constants.hpp>

using namespace libmpdataxx;

const int nt = 1400;

int main() 
{
//<listing-1>
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqns = 2 };
    enum { rhs_scheme = 
      solvers::rhs_scheme_t::trapez };
    struct ix { enum {psi, phi}; };
  };
//</listing-1>

  using real_t = typename ct_params_t::real_t;

  using sim_t = output::gnuplot<
    coupled_harmosc_stats<ct_params_t>
  >;
  typename sim_t::rt_params_t p; 

//<listing-2>
  // run-time parameters
  using boost::math::constants::pi;
  p.dt = 1;
  p.omega = 2 * pi<real_t>() / p.dt / 400;
//</listing-2>
  p.grid_size = {1001};
  p.outfreq = 10; 

  using ix = typename ct_params_t::ix;
  p.outvars = {
    {ix::psi, {.name = "psi", .unit = "1"}},
    {ix::phi, {.name = "phi", .unit = "1"}}
  };
  p.gnuplot_command = "plot";

  // instantiation
  concurr::threads<sim_t, bcond::cyclic, bcond::cyclic> run(p);

  // initial condition
  {
    blitz::firstIndex i;
    run.advectee(ix::psi) = where(        
      i<= 50 || i>= 150,                            // if
      0,                                            // then
      0.5 * (1 + cos(2 * pi<real_t>() * i / 100))   // else
    );
    run.advectee(ix::phi) = real_t(0);
  }
  run.advector() = .5;

  // integration
  run.advance(nt);
}
