/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "shallow_water.hpp"
#include <libmpdata++/concurr/serial.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx; 

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

const int nx = 40, nt = 0;

int main() 
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 1 };
    enum { n_eqs = 2 };
 
    // options
    enum { rhs_scheme = solvers::strang };

    // indices
    struct ix { enum {qx, h, vip_i=qx, vip_den=h}; };

    // hints
    enum { hint_norhs = formulae::opts::bit(ix::h) }; 
    // enum { hint_noneg = formulae::opts::bit(ix::h) };  // TODO: reconsider?
  };
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;

  // solver & output choice
  using solver_t = output::gnuplot<
    shallow_water<ct_params_t>
  >;

  // run-time parameters
  solver_t::rt_params_t p; 

  p.span = { nx };

  p.dt = .05;
  p.di = 1;

  //p.g = 1;
  p.outfreq = 1; //nt / 25;
  p.outvars = {
    {ix::h, { .name = "h",   .unit = "m" }}, 
  };
  p.gnuplot_command = "plot";
  p.gnuplot_with = "steps";

  // instantiation
  concurr::serial<solver_t, bcond::cyclic, bcond::cyclic> run(p);

  // initial condition
  const real_t dx = .1;
  {
    blitz::firstIndex i;

    struct intcond
    {
      real_t operator()(const real_t &x) const
      {
	 return 
	   std::abs(x) <= 1 // if
	   ? 1 - x*x        // then
	   : 0;             // else
      }
      BZ_DECLARE_FUNCTOR(intcond);
    };

    run.advectee(ix::h) = intcond()(dx*(i+.5) - nx*dx/2);

    run.advectee(ix::qx) = 0;
  }

  // integration
  run.advance(nt); 
};
