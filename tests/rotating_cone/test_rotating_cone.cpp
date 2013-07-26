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

#include <libmpdata++/solvers/adv/mpdata_2d.hpp>
#include <libmpdata++/solvers/adv/mpdata_fct_2d.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

enum {x, y};
using real_t = float;

real_t 
  dt = .1,
  dx = 1,
  dy = 1,
  omega = -.1,
  h = 4., // TODO: other name!
  h0 = 100.; // change it to 1 to see scary things!

// Anderson-Fattachi?
//    dt = 10 * pi<real_t>(),
//    omega = -.001,// / (2 * pi<real_t>()),
//    r = 4. * dx,
//    h0 = -.5,
//    x0 = 21. * dx,
//    y0 = 15. * dy;

/// @brief settings from @copybrief Anderson_and_Fattahi_1974
template <class T>
void setup(T &solver, int n[2]) 
{
  real_t
    r = 15. * dx,
    x0 = 75 * dx,
    y0 = 50 * dy,
    xc = .5 * n[x] * dx,
    yc = .5 * n[y] * dy;

  blitz::firstIndex i;
  blitz::secondIndex j;

  // cone shape
  decltype(solver.state()) tmp(solver.state().extent());
  tmp = blitz::pow((i+.5) * dx - x0, 2) + blitz::pow((j+.5) * dy - y0, 2);
  solver.state() = h0 + where(tmp - pow(r, 2) <= 0, h * blitz::sqr(1 - tmp / pow(r, 2)), 0.);

  // constant angular velocity rotational field
  solver.courant(x) = -omega * (j * dy - yc) * dt / dx;
  solver.courant(y) =  omega * (i * dx - xc) * dt / dy;
}

template <class T>
void setopts(T &p, int nt, int n_iters)
{
  p.outfreq = nt; //nt;///10; // TODO
  p.outvars[0].name = "psi";
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << n_iters << "_%s_%d.svg";
    p.gnuplot_output = tmp.str();    
  }
  p.gnuplot_view = "map";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
  p.gnuplot_contour = true;
  {
    std::ostringstream tmp;
    tmp << "[" << h0 -.5 << " : " << h0 + h + .5 << "]";
    p.gnuplot_cbrange = tmp.str();
  }
  p.gnuplot_xrange = "[60 : 90]";
  p.gnuplot_yrange = "[35 : 65]";
  p.gnuplot_maxcolors = 10;
  {
    std::ostringstream tmp;
    tmp << "levels incremental " << h0 -.25 << ", .25," << h0 + h + .25;
    p.gnuplot_cntrparam = tmp.str();
  }
  p.gnuplot_term = "svg";

/*
  p.outfreq = nt;
  p.gnuplot_with = "lines";
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << n_iters << "_%s_%d.svg";
    p.gnuplot_output = tmp.str();
  }
  p.outvars[0].name = "psi";
*/
}

int main() 
{
  using namespace libmpdataxx;

//  int n[] = {32, 32}, nt = 200;
  int n[] = {101, 101}, nt = 628 * 6;

  const int n_iters = 2;
  //using solver_t = output::gnuplot<solvers::mpdata_2d<real_t, n_iters, 1/*, formulae::mpdata::toa*/>>;
  using solver_t = output::gnuplot<solvers::mpdata_fct_2d<real_t, n_iters, 1, formulae::mpdata::iga | formulae::mpdata::toa>>;
  solver_t::params_t p;
  setopts(p, nt, n_iters);
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(n[x], n[y], p); 

  setup(slv, n); 
  slv.advance(nt);
}

