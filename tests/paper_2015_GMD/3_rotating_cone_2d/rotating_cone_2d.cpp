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
#include "rotating_cone_stats.hpp"

using namespace libmpdataxx;

template <int opts_arg, int opts_iters>
void test(const std::string filename)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
//<listing-1>
    enum { n_dims = 2 };
    enum { n_eqns = 1 };
//</listing-1>
    enum { opts = opts_arg };
  };

  typename ct_params_t::real_t    /// @brief settings from @copybrief Anderson_and_Fattahi_1974
    dt = .1,                      //    dt = 10 * pi<real_t>(),
    dx = 1,                       //    x0 = 21. * dx,
    dy = 1,                       //    y0 = 15. * dy;
    omega = .1,                   //    omega = -.001,// / (2 * pi<real_t>()),
    h = 4.,                       //    r = 4. * dx,
    h0 = 1;                       //    h0 = -.5,

  int nt = 628 * 6;

  using slv_out_t = 
    stats<
      output::gnuplot<
        solvers::mpdata<ct_params_t>
      >
    >;
  typename slv_out_t::rt_params_t p;
  p.dt = dt;  //to have this->dt in stats 

  // pre instantiation
  switch (opts_iters) // the crazy logic below is just for prettying the listing!
  {
    case 3: 
//<listing-4>
      p.n_iters = 3;
//</listing-4>
      break;
    default:
      p.n_iters = opts_iters; 
  }
  p.grid_size = {101, 101};

  p.outfreq = nt * dt; 
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
  {
    std::ostringstream tmp;
    tmp << "levels incremental " << h0 -.25 << ", .25," << h0 + h + .25;
    p.gnuplot_cntrparam = tmp.str();
  }
  p.gnuplot_fontsize = "14";
  p.gnuplot_cbrange = "[.75 : 5.25]";
  p.gnuplot_palette = "defined (" 
    "0.75 '#ff0000',"
    "1.00 '#ff0000',"
    "1.00 '#ffffff',"
    "1.25 '#ffffff',"
    "1.25 '#993399',"
    "2.25 '#00CCFF',"
    "3.25 '#66CC00',"
    "4.25 '#FC8727',"
    "5.25 '#FFFF00') maxcolors 18";
  p.gnuplot_term = "svg";
  p.gnuplot_title = "notitle";

//<listing-2>
  // instantiation
  concurr::threads<
    slv_out_t, 
    bcond::open, bcond::open,
    bcond::open, bcond::open
  > run(p); 
//</listing-2>
  {

    // constants used in the set-up definition
    enum {x, y};
    const typename ct_params_t::real_t
      r = 15. * dx,
      x0 = 50 * dx,
      y0 = 75 * dy,
      xc = .5 * (p.grid_size[x]-1) * dx,
      yc = .5 * (p.grid_size[y]-1) * dy;

//<listing-3>
    // temporary array of the same ...
    decltype(run.advectee())        // type 
      tmp(run.advectee().extent()); // and size 
    // ... as the one returned by advectee()

    // helper vars for Blitz++ tensor notation
    blitz::firstIndex i;
    blitz::secondIndex j;

    // cone shape ...
    tmp = blitz::pow(i * dx - x0, 2) + 
          blitz::pow(j * dy - y0, 2);

    // ... cut off at zero
    run.advectee() = h0 + where(
      tmp - pow(r, 2) <= 0,                // if
      h - blitz::sqrt(tmp / pow(r/h,2)),   // then  
      0.                                   // else
    );

    // constant-angular-velocity rotational field
    run.advector(x) =  omega * (j * dy - yc) * dt/dx;
    run.advector(y) = -omega * (i * dx - xc) * dt/dy;
//</listing-3>
  }
  // TODO: an assert confirming that the above did what it should have done
  //       (in context of the advector()'s use of blitz::Array::reindex())

  // time stepping
  run.advance(nt * dt);
}

int main()
{
  {
    enum { opts = 0 };
    enum { opts_iters = 2};
    test<opts, opts_iters>("basic");
  }
  {
    enum { opts = opts::fct };
    enum { opts_iters = 2};
    test<opts, opts_iters>("fct");
  }
  {
    enum { opts = opts::fct | opts::tot };
    enum { opts_iters = 3};
    test<opts, opts_iters>("iters3_tot_fct");
  }
  {
    enum { opts = opts::iga | opts::fct};
    enum { opts_iters = 2};
    test<opts, opts_iters>("iga_fct");
  }
  {
    enum { opts = opts::iga | opts::tot | opts::fct };
    enum { opts_iters = 2};
    test<opts, opts_iters>("iga_tot_fct");
  }
}
