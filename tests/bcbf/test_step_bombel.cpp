/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief test for pressure solvers (as in Smolarkiewicz & Pudykiewicz 1992, fig.3) 
 * buoyant convection in Boussinesq flow
 */

#include "bombel.hpp"
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx;

int main() 
{
  // compile-time parameters
//<listing-1>
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqns = 3 };
    enum { rhs_scheme = solvers::euler_b };
    enum { prs_scheme = solvers::cr };
    struct ix { enum {
      u, w, tht, 
      vip_i=u, vip_j=w, vip_den=-1
    }; };
  }; 
//</listing-1>
  using ix = typename ct_params_t::ix;

  const int r0 = 250; 
  const int nx = 200, ny = 200, nt = 800;
  typename ct_params_t::real_t Tht_amb = 1; //1; //300; // ambient state (constant thoughout the domain)

  // conjugate residual
  using solver_t = output::gnuplot<bombel<ct_params_t>>;

  // run-time parameters
  solver_t::rt_params_t p;

  p.dt = .75;
  p.di = p.dj = 10.; 
  p.Tht_amb = Tht_amb; 

  p.outfreq = 100; //12;
  p.outvars = {
//    {ix::u,   {.name = "u",   .unit = "m/s"}}, 
//    {ix::w,   {.name = "w",   .unit = "m/s"}}, 
    {ix::tht, {.name = "tht", .unit = "K"  }}
  };
  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_%s_%04d.svg";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
  p.gnuplot_contour = true;
//  rt_params.gnuplot_cbrange = "[299.85 : 300.65]";
  p.gnuplot_cbrange = "[299.85 - 299 : 300.65 - 299]";
  p.gnuplot_maxcolors = 8;
//  rt_params.gnuplot_cntrparam = "levels incremental 299.85, 0.1, 300.65";
  p.gnuplot_cntrparam = "levels incremental 299.85 - 299, 0.1, 300.65 - 299";
  p.gnuplot_term = "svg";
//<listing-2>
  p.tol = 1e-5;
//</listing-2>
  p.grid_size = {nx, ny};

  libmpdataxx::concurr::threads<
    solver_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > slv(p);

  {
    // initial condition
    blitz::firstIndex i;
    blitz::secondIndex j;

    slv.advectee(ix::tht) = Tht_amb + where(
      // if
      pow(i * p.di - 4    * r0 , 2) + 
      pow(j * p.dj - 1.04 * r0 , 2) <= pow(r0, 2), 
      // then
      .5, 
      // else
      0
    );
    slv.advectee(ix::u) = 0; 
    slv.advectee(ix::w) = 0; 
  }

  // integration
  slv.advance(nt); 
};
