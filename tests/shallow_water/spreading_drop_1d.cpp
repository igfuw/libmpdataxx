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
  enum { opts = opts_arg };
  enum { rhs_scheme = solvers::trapez };
  
  //enum { fp_round_mode = FE_TONEAREST };

  // indices
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

// TODO: all this plotting logic should be done with a new libmpdataxx::output
struct out_t
{
  std::ofstream x, h, q, t;
  out_t(const std::string &pfx) : 
    x(pfx + ".x"), 
    h(pfx + ".h"), 
    q(pfx + ".q"), 
    t(pfx + ".t")
  {}
};

template <class ix, class run_t>
void output(run_t &run, const int &t, const real_t &dx, const real_t &dt, out_t &out)
{
  // x coordinate (once)
  if (t == 0)
  {
    for (int i = 0; i < run.advectee().extent(0); ++i) 
      out.x << i * dx << "\t";
    out.x << "\n";
  } 

  // time
  out.t << t * dt << "\t" << "\n"; 
  
  // layer depth
  for (auto &it : run.advectee(ix::h)) out.h << it << "\t";
  out.h << "\n";

  // momentum
  for (auto &it : run.advectee(ix::qx)) out.q << it << "\t";
  out.q << "\n";
}

template<int opts>
void test(const std::string &pfx) 
{
  using ix = typename ct_params_t<opts>::ix;

  // solver choice
  using solver_t = shallow_water<ct_params_t<opts>>;

  // run-time parameters
  typename solver_t::rt_params_t p; 

  p.dt = .01;
  p.di = .05;
  p.grid_size = { int(16 / p.di) };
  p.g = 1;
  p.vip_eps = 1.e-10; // in 1D apparently it's enough!

  // instantiation
  concurr::serial<
    solver_t, 
    bcond::cyclic, bcond::cyclic
  > run(p); // TODO: change into open bc

  // initial condition
  {
    blitz::firstIndex i;
    run.advectee(ix::h) = intcond()(p.di * i - (p.grid_size[0]-1) * p.di / 2);
  }
  run.advectee(ix::qx) = 0;

  // integration
  out_t out(pfx);
  output<ix>(run, 0, p.di, p.dt, out);
  for (int t = 0; t < nt; t += outfreq)
  {
    run.advance(outfreq); 
    output<ix>(run, t + outfreq, p.di, p.dt, out);
  }
}

int main()
{
  test<formulae::opts::iga | formulae::opts::fct>("fct+iga");
  test<formulae::opts::abs | formulae::opts::fct>("fct+abs");
  //plotting model results and analitic solution; 
  //python uses sys.argv[1:0] for choosing model outputs
  system("python ../../../tests/shallow_water/papierplot_shallow_water_1d.py fct+abs fct+iga ");
}

