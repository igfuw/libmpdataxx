/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * \include "rotating_cone/test_rotating_cone.cpp"
 * \image html "../../tests/rotating_cone/figure.svg"
 */

#include <cmath>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <libmpdata++/solvers/adv/mpdata_fct_2d.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

enum {x, y};
using real_t = double;

/// @brief settings from @copybrief Anderson_and_Fattahi_1974
template <class T>
void setup(T &solver, int n[2]) 
{
  real_t 
    dt = 10 * pi<real_t>(),
    dx = 1,
    dy = 1,
    xc = .5 * n[x] * dx,
    yc = .5 * n[y] * dy,
    omega = -.001,// / (2 * pi<real_t>()),
    h = 1., // TODO: other name!
    r = 4. * dx,
    h0 = -.5,
    x0 = 21. * dx,
    y0 = 15. * dy;

  blitz::firstIndex i;
  blitz::secondIndex j;

  // cone shape
  decltype(solver.state()) tmp(solver.state().extent());
  tmp = h - blitz::sqr(
    (
      blitz::pow((i+.5) * dx - x0, 2) +  
      blitz::pow((j+.5) * dy - y0, 2) 
    ) /
    std::pow(r / h, 2)
  );
  solver.state() = h0 + where(tmp > 0, tmp, 0.);

  // constant angular velocity rotational field
  solver.courant(x) = -omega * (j * dy - yc) * dt / dx;
  solver.courant(y) = omega * (i * dx - xc) * dt / dy;
}

template <class T>
void setopts(T &p, int nt, int n_iters)
{
  //p.outfreq = nt/10; // TODO
  p.gnuplot_with = "lines";
  p.gnuplot_zrange = p.gnuplot_cbrange = "[-1:1]";
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << n_iters << "_%s_%d.svg";
    p.gnuplot_output = tmp.str();    
  }
  p.outvars = {{0, {.name = "psi", .unit = "1"}}};
}

int main() 
{
  using namespace libmpdataxx;

  int n[] = {32, 32}, nt = 200;

  const int n_iters = 2;
  using solver_t = output::gnuplot<solvers::mpdata_fct_2d<real_t, n_iters>>;
  solver_t::params_t p;
  setopts(p, nt, n_iters);
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p); 

  setup(slv, n); 
  slv.advance(nt);
}
