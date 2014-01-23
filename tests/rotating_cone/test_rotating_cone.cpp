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

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx;

enum {x, y};
struct ct_params_t : ct_params_default_t
{
  using real_t = float;
  enum { n_dims = 2 };
  enum { n_eqs = 1 };
  enum { opts = formulae::opts::iga | formulae::opts::tot | formulae::opts::fct };
};

ct_params_t::real_t 
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

int main() 
{
  int nt = 628 * 6;

  using solver_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  solver_t::rt_params_t p;

  // pre instantiation
  p.n_iters = 2;
  p.span = {101, 101};

  p.outfreq = nt; 
  p.outvars[0].name = "psi";
  {
    std::ostringstream tmp;
    tmp << "figure_iters=" << p.n_iters << "_%s_%d.svg";
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

  // instantiation
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic> slv(p); 

  // post instantiation
  {
    ct_params_t::real_t
      r = 15. * dx,
      x0 = 75 * dx,
      y0 = 50 * dy,
      xc = .5 * p.span[x] * dx,
      yc = .5 * p.span[y] * dy;

    blitz::firstIndex i;
    blitz::secondIndex j;

    // cone shape
    decltype(slv.advectee()) tmp(slv.advectee().extent());
    tmp = blitz::pow((i+.5) * dx - x0, 2) + blitz::pow((j+.5) * dy - y0, 2);
    slv.advectee() = h0 + where(tmp - pow(r, 2) <= 0, h * blitz::sqr(1 - tmp / pow(r, 2)), 0.);

    // constant angular velocity rotational field
    slv.advector(x) = -omega * (j * dy - yc) * dt / dx;
    slv.advector(y) =  omega * (i * dx - xc) * dt / dy;
    // TODO: an assert confirming that the above did what it should have done
    //       (in context of the advector()'s use of blitz::Array::reindex())
  }

  // time stepping
  slv.advance(nt);
}

