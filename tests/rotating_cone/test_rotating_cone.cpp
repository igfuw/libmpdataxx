/* 
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

template <int opts_arg, int opts_iters>
void test(const std::string filename)
{
  enum {x, y};
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
//<listing-1>
    enum { n_dims = 2 };
    enum { n_eqns = 1 };
//</listing-1>
    enum { opts = opts_arg };
  };

  typename ct_params_t::real_t 
    dt = .1,
    dx = 1,
    dy = 1,
    omg = .1,
    h = 4., // TODO: other name!
    h0 = 1;
    //  h0 = 100.; // change it to 1 to see scary things!

/// @brief settings from @copybrief Anderson_and_Fattahi_1974
//    dt = 10 * pi<real_t>(),
//    omg = -.001,// / (2 * pi<real_t>()),
//    r = 4. * dx,
//    h0 = -.5,
//    x0 = 21. * dx,
//    y0 = 15. * dy;

  int nt = 628 * 6;

  using sim_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  typename sim_t::rt_params_t p;

  // pre instantiation
  p.n_iters = opts_iters; 
  p.grid_size = {101, 101};

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
//  p.gnuplot_xrange = "[0 : 100]";
//  p.gnuplot_yrange = "[0 : 100]";
  p.gnuplot_maxcolors = 10;
  {
    std::ostringstream tmp;
    tmp << "levels incremental " << h0 -.25 << ", .25," << h0 + h + .25;
    p.gnuplot_cntrparam = tmp.str();
  }
  p.gnuplot_term = "svg";

//<listing-2>
  // instantiation
  concurr::threads<
    sim_t, 
    bcond::open, bcond::open,
    bcond::open, bcond::open
  > run(p); 
//</listing-2>
  {
//<listing-3>
    // post instantiation
    typename ct_params_t::real_t
      r = 15. * dx,
      x0 = 50 * dx,
      y0 = 75 * dy,
      xc = .5 * (p.grid_size[x]-1) * dx,
      yc = .5 * (p.grid_size[y]-1) * dy;

    blitz::firstIndex i;
    blitz::secondIndex j;

    // cone shape
    decltype(run.advectee()) 
      tmp(run.advectee().extent());

    tmp = blitz::pow(i * dx - x0, 2) + 
          blitz::pow(j * dy - y0, 2);

    run.advectee() = h0 + where(
      tmp - pow(r, 2) <= 0,                  //if
      h * blitz::sqr(1 - tmp / pow(r, 2)),   //then
      0.                                     //else
    );

    // constant angular velocity rotational field
    run.advector(x) =  omg * (j * dy - yc) * dt/dx;
    run.advector(y) = -omg * (i * dx - xc) * dt/dy;
//</listing-3>
  }
    // TODO: an assert confirming that the above did what it should have done
    //       (in context of the advector()'s use of blitz::Array::reindex())

  // time stepping
  run.advance(nt);
  
  std::cout<<"min(psi) = " << min(run.advectee()) << std::endl;
}

int main()
{
  {
    enum { opts = 0 };
    const int opts_iters = 2;
    test<opts, opts_iters>("basic");
  }
  {
    enum { opts = formulae::opts::fct };
    const int opts_iters = 2;
    test<opts, opts_iters>("fct");
  }
  {
    enum { opts = formulae::opts::fct | formulae::opts::tot };
    const int opts_iters = 3;
    test<opts, opts_iters>("iters3_tot_fct");
  }
  {
    enum { opts = formulae::opts::iga | formulae::opts::fct};
    const int opts_iters = 2;
    test<opts, opts_iters>("iga_fct");
  }
  {
    enum { opts = formulae::opts::iga | formulae::opts::tot | formulae::opts::fct };
    const int opts_iters = 2;
    test<opts, opts_iters>("iga_tot_fct");
  }
}
