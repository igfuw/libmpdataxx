/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#include "shallow_water.hpp"
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/hdf5_xdmf.hpp>
using namespace libmpdataxx; 

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <fstream>

const int 
  nt = 300,
  outfreq = 100;

// compile-time parameters
// enum { hint_noneg = formulae::opts::bit(ix::h) };  // TODO: reconsider?
//<listing-1>
struct ct_params_t : ct_params_default_t
{
  using real_t = double;
  enum { n_dims = 2 };
  enum { n_eqs = 3 };
  
  // options
  enum { opts = formulae::opts::fct | formulae::opts::iga };
  enum { rhs_scheme = solvers::strang };
  
  // indices
  struct ix { enum {
      qx, qy, h, 
      vip_i=qx, vip_j=qy, vip_den=h
  }; }; 
  
  // hints
  enum { hint_norhs = formulae::opts::bit(ix::h) }; 
};
//</listing-1>

using real_t = typename ct_params_t::real_t;
using ix = typename ct_params_t::ix;

struct intcond
{
  real_t operator()(const real_t &x, const real_t &y) const
  {
    return 
      x*x + y*y <= 1 // if
      ? 1 - x*x - y*y  // then
      : 0;             // else
  }
  BZ_DECLARE_FUNCTOR2(intcond);
};


int main() 
{
  using ix = typename ct_params_t::ix;

  // solver choice
  using solver_t = output::hdf5_xdmf<shallow_water<ct_params_t>>;


  // run-time parameters
  solver_t::rt_params_t p; 

  p.dt = .01;
  p.di = .05;
  p.dj = .05;
  p.span = { int(16 / p.di), int(16 / p.dj) };
  p.g = 1;
  p.outfreq = outfreq;
  p.outdir = "spreading_drop_2d.out";
  p.outvars = {
    {ix::h,  {.name="h" , .unit="TODO"}}, 
    {ix::qx, {.name="qx", .unit="TODO"}}, 
    {ix::qy, {.name="qy", .unit="TODO"}}
  };
  //p.vip_eps = 1e-5;

  // instantiation
  concurr::threads<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > run(p); // TODO: change into open bc

  // initial condition
  {
    blitz::firstIndex i;
    blitz::secondIndex j;
    run.advectee(ix::h) = intcond()(
      p.di * i - (p.span[0]-1) * p.di / 2, 
      p.dj * j - (p.span[1]-1) * p.dj / 2
    );
  }
  run.advectee(ix::qx) = 0;
  run.advectee(ix::qy) = 0;

  run.advance(nt);
};
