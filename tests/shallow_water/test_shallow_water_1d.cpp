/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "shallow_water.hpp"
#include <libmpdata++/concurr/serial.hpp>
using namespace libmpdataxx; 

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <fstream>

const int 
  nt = 300,
  outfreq = 1;

// compile-time parameters
// enum { hint_noneg = formulae::opts::bit(ix::h) };  // TODO: reconsider?
//<listing-1>
struct ct_params_t : ct_params_default_t
{
  using real_t = double;
  enum { n_dims = 1 };
  enum { n_eqs = 2 };
  
  // options
  enum { opts = 0 };
  enum { rhs_scheme = solvers::strang };
  
  // indices
  struct ix { enum {
    qx, h, 
    vip_i=qx, vip_den=h
  }; }; 
  
  // hints
  enum { hint_norhs = formulae::opts::bit(ix::h) }; 
};
//</listing-1>

using real_t = typename ct_params_t::real_t;
using ix = typename ct_params_t::ix;

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

std::ofstream ox("out.x"), oh("out.h"), oq("out.q"), ot("out.t");

template <class run_t>
void output(run_t &run, const int &t, const real_t &dx, const real_t &dt)
{
  // x coordinate (once)
  if (t == 0)
  {
    for (int i = 0; i < run.advectee().extent(0); ++i) 
      ox << (i + .5) * dx << "\t";
    ox << "\n";
  } 

  // time
  ot << t * dt << "\t" << "\n"; 
  
  // layer depth
  for (auto &it : run.advectee(ix::h)) oh << it << "\t";
  oh << "\n";

  // momentum
  for (auto &it : run.advectee(ix::qx)) oq << it << "\t";
  oq << "\n";
}

int main() 
{
  using ix = typename ct_params_t::ix;

  // solver choice
  using solver_t = shallow_water<ct_params_t>;

  // run-time parameters
  solver_t::rt_params_t p; 

  p.dt = .01;
  p.di = .05;
  p.span = { int(16 / p.di) };
  p.g = 1;

  // instantiation
  concurr::serial<solver_t, bcond::cyclic, bcond::cyclic> run(p);

  // initial condition
  {
    blitz::firstIndex i;
    run.advectee(ix::h) = intcond()(p.di * (i+.5) - p.span[0] * p.di / 2);
  }
  run.advectee(ix::qx) = 0;

  // integration
  output(run, 0, p.di, p.dt);
  for (int t = 0; t < nt; t += outfreq)
  {
    run.advance(outfreq); 
    output(run, t + outfreq, p.di, p.dt);
  }
};
