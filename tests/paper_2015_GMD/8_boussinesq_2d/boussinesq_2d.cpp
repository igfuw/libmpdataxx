/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief test for pressure solvers (as in Smolarkiewicz & Pudykiewicz 1992, fig.3) 
 * buoyant convection in Boussinesq flow
 */
#include <libmpdata++/solvers/boussinesq.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
#include "boussinesq_stats.hpp"

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
  using real_t = typename ct_params_t::real_t;

  const int r0 = 250; 
  const int nx = 201, ny = 201, nt = 800;
  typename ct_params_t::real_t Tht_ref = 300; // reference state (constant throughout the domain)

  // conjugate residual
  using slv_out_t = 
    stats< 
      output::gnuplot<
        solvers::boussinesq<ct_params_t>
      >
    >;
  // run-time parameters
  slv_out_t::rt_params_t p;

  p.dt = .75;
  p.di = p.dj = 10.; 
  p.Tht_ref = Tht_ref; 

  p.outfreq = 100; //12;
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

  real_t eps = .01;

  p.gnuplot_cntrparam = "levels incremental 299.95, 0.1, 300.55";
  p.gnuplot_cbrange = "[299.95 : 300.55]";
  p.gnuplot_cbtics = "300.05, 0.1, 300.45";
  p.gnuplot_palette = "defined ("
    "299.95 '#ffffff', "
    "300.05 '#ffffff', 300.05 '#993399', "
    "300.15 '#993399', 300.15 '#00CCFF', "
    "300.25 '#00CCFF', 300.25 '#66CC00', "
    "300.35 '#66CC00', 300.35 '#FC8727', "
    "300.45 '#FC8727', 300.45 '#FFFF00', "
    "300.55 '#FFFF00'"
  ")";
  p.gnuplot_term = "svg";
  p.prs_tol = 1e-7;
  p.grid_size = {nx, ny};

  libmpdataxx::concurr::threads<
    slv_out_t, 
    bcond::cyclic, bcond::cyclic,
    bcond::cyclic, bcond::cyclic
  > slv(p);

  {
    // initial condition
    blitz::firstIndex i;
    blitz::secondIndex j;
    
    slv.sclr_array("tht_e") = Tht_ref;
    slv.advectee(ix::tht) = Tht_ref + where(
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

  if (min(slv.advectee(ix::tht)) < 300-eps || max(slv.advectee(ix::tht)) > 300.5+eps)
    throw std::runtime_error("too big under- or over-shots :("); 
};
