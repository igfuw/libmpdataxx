/* 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
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
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqns = 1 };
    enum { opts = opts_arg };
  };

  typename ct_params_t::real_t 
    dt = .1,
    dx = 1,
    dy = 1,
    omg = .1,
    h = 4., // TODO: other name!
    h0 = 1, 
    h0d = 1,
    hd = 10;
  
int nt = 628 * 2;

  using slv_out_t = output::gnuplot<solvers::mpdata<ct_params_t>>;
  typename slv_out_t::rt_params_t p;

  // pre instantiation
  p.n_iters = opts_iters; 
  p.grid_size = {101, 101};

  p.outfreq = 10; 
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
//  p.gnuplot_xrange = "[25 : 75]";
//  p.gnuplot_yrange = "[50 : 100]";
  p.gnuplot_xrange = "[0 : 100]";
  p.gnuplot_yrange = "[0 : 100]";
  p.gnuplot_maxcolors = 10;
  {
    std::ostringstream tmp;
    tmp << "levels incremental " << h0 -.25 << ", .25," << h0 + h + .25;
    p.gnuplot_cntrparam = tmp.str();
  }
  p.gnuplot_fontsize = "14";
  p.gnuplot_term = "svg";
  p.gnuplot_title = "notitle";

  // instantiation
  concurr::threads<
    slv_out_t, 
    bcond::open, bcond::open,
    bcond::open, bcond::open
  > run(p); 
  {

    // constants used in the set-up definition
    enum {x, y};
    const typename ct_params_t::real_t
      r = 15. * dx,
      rd = 10 * dx,
      x0 = 50 * dx,
      y0 = 75 * dy,
      x0d = 25 * dx,    //location of desity perturbation
      y0d = 50 * dy,
      xc = .5 * (p.grid_size[x]-1) * dx,
      yc = .5 * (p.grid_size[y]-1) * dy;

    // temporary array of the same ...
    decltype(run.advectee())        // type 
      tmp(run.advectee().extent()), // and size
      density_tmp(run.advectee().extent()); 
    // ... as the one returned by advectee()

    // helper vars for Blitz++ tensor notation
    blitz::firstIndex i;
    blitz::secondIndex j;

    // cone shape ...
    tmp = blitz::pow(i * dx - x0, 2) + 
          blitz::pow(j * dy - y0, 2);

    // ... cut off at zero
    run.advectee() = h0 + where(
      tmp - pow(r, 2) <= 0,                  //if
      h * blitz::sqr(1 - tmp / pow(r, 2)),   //then
      0.                                     //else
    );

    // density shape (cone again ...)
    tmp = blitz::pow(i * dx - x0d, 2) + 
          blitz::pow(j * dy - y0d, 2);
    
/*    density_tmp = h0d + where(
      tmp - pow(rd, 2) <= 0,                   //if
      hd * blitz::sqr(1 - tmp / pow(rd, 2)),   //then
      0.                                       //else
    );
*/

   density_tmp = 1.;

    // constant-angular-velocity rotational field
    run.advector(x) =  1. / density_tmp * omg * (j * dy - yc) * dt/dx;
    run.advector(y) = -1. / density_tmp * omg * (i * dx - xc) * dt/dy;
   
    run.g_factor() = density_tmp;

  }
  // time stepping
  run.advance(nt);
  
  std::cout<<"min(psi) = " << min(run.advectee()) << std::endl;
}

int main()
{
/*  {
    enum { opts = 0 | opts::nug };
    enum { opts_iters = 2};
    test<opts, opts_iters>("basic");
  }
*/
  {
    enum { opts = opts::fct | opts::nug };
    enum { opts_iters = 2};
    test<opts, opts_iters>("fct");
  }
/*
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
*/
}
