/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief test of fixed bcond
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
  using real_t = typename ct_params_t::real_t;

  const int nx = 401, ny = 101, nt = 2000;
  typename ct_params_t::real_t Tht_ref = 300; // reference state (constant throughout the domain)

  // conjugate residual
  using slv_out_t = 
    output::gnuplot<
      solvers::boussinesq<ct_params_t>
  >;
  // run-time parameters
  slv_out_t::rt_params_t p;

  p.dt = .75;
  p.di = p.dj = 10.; 
  p.Tht_ref = Tht_ref; 

  p.outfreq = 200; //12;
  p.outvars = {
   {ix::u,   {.name = "u",   .unit = "m/s"}}, 
   {ix::w,   {.name = "w",   .unit = "m/s"}}, 
    {ix::tht, {"tht", "K"  }}
  };
  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_%s_%04d.svg";
  p.gnuplot_size = "ratio " + std::to_string(double(ny) / nx);

  real_t eps = .01;

  const double term_x_size = 1200; // in pixels
  const double term_y_size = term_x_size * double(ny) / nx;
  p.gnuplot_term = "svg size " + std::to_string(int(term_x_size)) + "," + std::to_string(int(term_y_size));
  p.prs_tol = 1e-7;
  p.grid_size = {nx, ny};

  libmpdataxx::concurr::threads<
    slv_out_t, 
    bcond::fixed, bcond::open,
    bcond::open, bcond::open
  > slv(p);

  {
    // initial condition
    blitz::firstIndex i;
    blitz::secondIndex j;
    
    slv.sclr_array("tht_e") = Tht_ref;
    slv.advectee(ix::tht) = Tht_ref
     + where(
      // if
      i < 1 && j > 11 && j < 15, 
      // then
      1, 
      // else
      0
    );
    slv.advectee(ix::u) = 2; 
    slv.advectee(ix::w) = 0; 
  }

  // integration
  slv.advance(nt);  
};
