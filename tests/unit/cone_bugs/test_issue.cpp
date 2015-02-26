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

#include <libmpdata++/solvers/mpdata.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

/// @brief settings from @copybrief Anderson_and_Fattahi_1974

using namespace libmpdataxx;

enum {x, y};

struct ct_params_t :ct_params_default_t
{
  using real_t = float;
  enum { n_dims = 2 };
  enum { n_eqns = 1 };
  enum { opts = opts::fct };
};

ct_params_t::real_t
  dt = .1,
  dx = 1,
  dy = 1,
  h = 4., // TODO: other name!
  h0 = 1.;

int main() 
{
  const int nt = 2; 

  using solver_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  solver_t::rt_params_t p;

  // pre instantiatio
  {
    p.n_iters = 2;
    p.grid_size = {15, 14};

    p.outfreq = nt;
    p.outvars[0].name = "psi";
    {
      std::ostringstream tmp;
      tmp << "bug_iters=" << p.n_iters << "_%s_%d.svg";
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
    {
      std::ostringstream tmp;
      tmp << "levels incremental " << h0 -.25 << ", .25," << h0 + h + .25;
      p.gnuplot_cntrparam = tmp.str();
    }
    p.gnuplot_term = "svg";
  }

  // instantiation
  concurr::threads<solver_t, bcond::cyclic, bcond::cyclic, bcond::cyclic, bcond::cyclic> slv(p); 

  // post instantiation
  {
    ct_params_t::real_t
      r = 12. * dx,
      x0 = 1,//75 * dx,
      y0 = 0, //50 * dy,
      xc = .5 * (p.grid_size[x]-1) * dx,
      yc = .5 * (p.grid_size[y]-1) * dy;

    blitz::firstIndex i;
    blitz::secondIndex j;

    // cone shape
    decltype(slv.advectee()) tmp(slv.advectee().extent());
    tmp = blitz::pow(i * dx - x0, 2) + blitz::pow(j * dy - y0, 2);
    slv.advectee() = h0 + where(tmp - pow(r, 2) <= 0, h * blitz::sqr(1 - tmp / pow(r, 2)), 0.);

    // constant angular velocity rotational field
    slv.advector(x) = .3; 
    slv.advector(y) = .3; 
  }

  // 
  slv.advance(nt);
  if (blitz::min(slv.advectee()) - 1 != 0) 
    throw std::runtime_error("values below one!");
}
