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

const int nt = 50;

// compile-time parameters
struct ct_params_t : ct_params_default_t
{
  using real_t = double;
  enum { n_dims = 1 };
  enum { n_eqs = 2 };
  
  // options
  enum { opts = formulae::opts::fct };
  enum { rhs_scheme = solvers::strang };
  
  // indices
  struct ix { enum {qx, h, vip_i=qx, vip_den=h}; };
  
  // hints
  enum { hint_norhs = formulae::opts::bit(ix::h) }; 
  // enum { hint_noneg = formulae::opts::bit(ix::h) };  // TODO: reconsider?
};

using real_t = typename ct_params_t::real_t;

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

int main() 
{
  using ix = typename ct_params_t::ix;

  // solver & output choice
  using solver_t = output::gnuplot<
    shallow_water<ct_params_t>
  >;

  // run-time parameters
  solver_t::rt_params_t p; 

  p.dt = .01;
  p.di = .05;
  p.span = { int(16 / p.di) };

  //p.g = 1;
  p.outfreq = nt / 3;
  p.outvars = {
    {ix::h, { .name = "h",   .unit = "m" }}, 
  };
  p.gnuplot_command = "plot";
  p.gnuplot_with = "steps";

  // instantiation
  concurr::serial<solver_t, bcond::cyclic, bcond::cyclic> run(p);

  // initial condition
  {
    blitz::firstIndex i;
  
    run.advectee(ix::h) = 1e-3 + intcond()(p.di * (i+.5) - p.span[0] * p.di / 2);

    run.advectee(ix::qx) = 0;
  }

  // integration
  run.advance(nt); 
};
