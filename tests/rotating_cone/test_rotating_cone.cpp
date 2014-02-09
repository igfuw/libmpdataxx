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

template <int opts_arg>
void test(const std::string filename)
{
  enum {x, y};
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqs = 1 };
    enum { opts = opts_arg };
  };

  typename ct_params_t::real_t 
    dt = .1,
    dx = 1,
    dy = 1,
    omega = .1,
    h = 4., // TODO: other name!
    h0 = 0;
    //  h0 = 100.; // change it to 1 to see scary things!

/// @brief settings from @copybrief Anderson_and_Fattahi_1974
//    dt = 10 * pi<real_t>(),
//    omega = -.001,// / (2 * pi<real_t>()),
//    r = 4. * dx,
//    h0 = -.5,
//    x0 = 21. * dx,
//    y0 = 15. * dy;

  int nt = 628 * 6;

  using sim_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  typename sim_t::rt_params_t p;

  // pre instantiation
  p.n_iters = 2; 
  p.span = {101, 101};

  p.outfreq = nt; 
  p.outvars[0].name = "psi";
  {
    std::ostringstream tmp;
    tmp << filename << "_%s_%d.svg";
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
  p.gnuplot_xrange = "[25 : 75]";
  p.gnuplot_yrange = "[50 : 100]";
  p.gnuplot_maxcolors = 10;
  {
    std::ostringstream tmp;
    tmp << "levels incremental " << h0 -.25 << ", .25," << h0 + h + .25;
    p.gnuplot_cntrparam = tmp.str();
  }
  p.gnuplot_term = "svg";

  // instantiation
  concurr::threads<sim_t, bcond::cyclic, bcond::cyclic> run(p); 

  // post instantiation
  {
    typename ct_params_t::real_t
      r = 15. * dx,
      x0 = 50 * dx,//75 in the article?
      y0 = 75 * dy,//50 in the article?
      xc = .5 * p.span[x] * dx,
      yc = .5 * p.span[y] * dy;

    blitz::firstIndex i;
    blitz::secondIndex j;

    // cone shape
    decltype(run.advectee()) tmp(run.advectee().extent());
    tmp = blitz::pow((i+.5) * dx - x0, 2) + blitz::pow((j+.5) * dy - y0, 2);
    run.advectee() = h0 + where(tmp - pow(r, 2) <= 0, h * blitz::sqr(1 - tmp / pow(r, 2)), 0.);

    // constant angular velocity rotational field
    run.advector(x) = -omega * ((j+.5) * dy - yc) * dt / dx;
    run.advector(y) =  omega * ((i+.5) * dx - xc) * dt / dy;
    // TODO: an assert confirming that the above did what it should have done
    //       (in context of the advector()'s use of blitz::Array::reindex())
  }

  // time stepping
  run.advance(nt);
}

int main()
{
  {
    enum { opts = 0 };
    test<opts>("basic");
  }
  {
    enum { opts = formulae::opts::fct };
    test<opts>("fct");
  }
  {
    enum { opts = formulae::opts::iga | formulae::opts::tot};
    test<opts>("iga_tot");
  }
  {
    enum { opts = formulae::opts::iga | formulae::opts::fct};
    test<opts>("iga_fct");
  }
  {
    enum { opts = formulae::opts::iga | formulae::opts::tot | formulae::opts::fct};
    test<opts>("iga_fct_tot");
  }
}
