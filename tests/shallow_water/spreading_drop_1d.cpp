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

using real_t = double;

// compile-time parameters
// enum { hint_noneg = formulae::opts::bit(ix::h) };  // TODO: reconsider?
//<listing-1>
template <int opts_arg>
struct ct_params_t : ct_params_default_t
{
  using real_t = ::real_t;
  enum { n_dims = 1 };
  enum { n_eqs = 2 };
  
  // options
  enum { opts = opts_arg};
  enum { rhs_scheme = solvers::trapez };
  
  // indices - TODO move vip to separate enum
  struct ix { enum {
    qx, h, 
    vip_i=qx, vip_den=h
  }; }; 
  
  // hints
  enum { hint_norhs = formulae::opts::bit(ix::h) }; 
};
//</listing-1>

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

template <class ix, class run_t>
void output(run_t &run, const int &t, const real_t &dx, const real_t &dt)
{
  // x coordinate (once)
  if (t == 0)
  {
    for (int i = 0; i < run.advectee().extent(0); ++i) 
      ox << i * dx << "\t";
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

template<int opts>
void test() 
{
  using ix = typename ct_params_t<opts>::ix;

  // solver choice
  using solver_t = shallow_water<ct_params_t<opts>>;

  // run-time parameters
  typename solver_t::rt_params_t p; 

  p.dt = .01;
  p.di = .05;
  p.span = { int(16 / p.di) };
  p.g = 1;
  p.vip_eps = 0; // in 1D apparently it's enough!

  // instantiation
  concurr::serial<
    solver_t, 
    bcond::cyclic, bcond::cyclic
  > run(p); // TODO: change into open bc

  // initial condition
  {
    blitz::firstIndex i;
    run.advectee(ix::h) = intcond()(p.di * i - (p.span[0]-1) * p.di / 2);
  }
  run.advectee(ix::qx) = 0;

  // integration
  output<ix>(run, 0, p.di, p.dt);
  for (int t = 0; t < nt; t += outfreq)
  {
    run.advance(outfreq); 
    output<ix>(run, t + outfreq, p.di, p.dt);
  }
}

int main()
{
  test<formulae::opts::iga | formulae::opts::fct>();
}

