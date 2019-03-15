/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 */

#include <libmpdata++/solvers/boussinesq.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx;

int main() 
{
  // compile-time parameters
  struct ct_params_t : ct_params_default_t
  {
    using real_t = double;
    enum { n_dims = 2 };
    enum { n_eqns = 3 };
    enum { rhs_scheme = solvers::trapez };
    enum { prs_scheme = solvers::cr };
    struct ix { enum {
      u, w, tht, 
      vip_i=u, vip_j=w, vip_den=-1
    }; };
  }; 

  using ix = typename ct_params_t::ix;

  const int r0 = 500; 
  const int nx = 201, ny = 201, nt = 1200;
  typename ct_params_t::real_t Tht_ref = 300; //1; // reference state (constant throughout the domain)

  // conjugate residual
  using solver_t = output::gnuplot<solvers::boussinesq<ct_params_t>>;

  // run-time parameters
  solver_t::rt_params_t p;

  p.dt = .75;
  p.di = p.dj = 10.; 
  p.Tht_ref = Tht_ref; 

  p.outfreq = 100;
  p.outvars = {
//    {ix::u,   {.name = "u",   .unit = "m/s"}}, 
//    {ix::w,   {.name = "w",   .unit = "m/s"}}, 
    {ix::tht, {"tht", "K"  }}
  };
  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_%s_%04d.svg";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
  p.gnuplot_contour = true;
  p.gnuplot_cntrparam = "levels incremental 299.95, 0.1, 300.65";
  p.gnuplot_cbrange = "[299.95 : 300.65]";
  p.gnuplot_cbtics = "('299.99' 299.99, '300.10' 300.1, '300.20' 300.2, '300.30' 300.3, '300.40' 300.4, '300.50' 300.5, '300.60' 300.6)";
  p.gnuplot_palette = "defined (" 
    "299.95 '#ff0000'," //         
    "299.99 '#ff0000'," // 
    "299.99 '#ffffff'," //         /\-
    "300.00 '#ffffff'," //        /  \-
    "300.00 '#ffffff'," //  -----/    \---
    "300.05 '#ffffff'," // -----/      \---___
    "300.05 '#993399'," //     /        \-     ---
    "300.20 '#00CCFF'," //    /          \-       ---
    "300.35 '#66CC00'," //   /____________\-
    "300.50 '#FC8727'," //
    "300.65 '#FFFF00') maxcolors 14";
  p.gnuplot_term = "svg";
  p.prs_tol = 1e-7;
  p.grid_size = {nx, ny};

  libmpdataxx::concurr::threads<
    solver_t, 
    bcond::cyclic , bcond::cyclic,
    bcond::rigid , bcond::rigid
  > slv(p);

  {
    // initial condition
    blitz::firstIndex i;
    blitz::secondIndex j;

    slv.advectee(ix::tht) = Tht_ref + where(
      // if
      pow(i * p.di - 1.04    * r0 , 2) + 
      pow(j * p.dj - 1.04 * r0 , 2) <= pow(r0, 2), 
      // then
      -.5, 
      // else
      0
    );
    slv.advectee(ix::u) = 0; 
    slv.advectee(ix::w) = 0; 
  }

  // integration
  slv.advance(nt); 
};
