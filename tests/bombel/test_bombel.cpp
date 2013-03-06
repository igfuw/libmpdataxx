/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a bombel 
 *
 * @section FIGURE
 *
 * \image html "../../tests/bombel/figure.svg"
 *
 */

// advection (<> should be used instead of "" in normal usage) 
#include "advoocat/solvers/mpdata_2d.hpp"
#include "advoocat/solvers/solver_pressure_mr.hpp"
#include "advoocat/solvers/solver_pressure_cr.hpp"
#include "advoocat/solvers/solver_pressure_pc.hpp"
#include "advoocat/bcond/cyclic_2d.hpp"
#include "advoocat/concurr/threads.hpp"
//theta->pressure
#include "advoocat/formulae/diagnose_formulae.hpp"
//physical constants
#include "advoocat/formulae/phc.hpp"

// plotting
#include "advoocat/output/gnuplot.hpp"

// auto-deallocating containers
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

enum {u, w, tht};  //eqations
enum {x, z};       //dimensions

using real_t = double;
const int n_iters = 2, n_eqs = 3;

using namespace advoocat;

#include "bombel.hpp"

int main() 
{
  const int nx = 500, ny = 500, nt = 5, n_out=1;
//  const int nx = 20, ny = 20, nt = 1, n_out=1;

  rng_t i(0, nx-1);
  rng_t j(0, ny-1);
  const real_t halo = 1;  

  // TODO p.Prs_amb = formulae::diagnose::p(p.Tht_amb);
  real_t dt = .1, dx = 1., dz = 1.; // timestep
  real_t Tht_amb = 300; // ambient state (constant thoughout the domain)

  boost::ptr_vector<concurr::any<real_t, 2>> slvs;

  { // minimum residual
    using solver_t = output::gnuplot<
      bombel<
	solvers::pressure_mr<
	  solvers::mpdata_2d<real_t, n_iters, n_eqs>,
	  u, w
        >
      >
    >;
    solver_t::params_t p; 

    p.dt = dt; 
    p.dx = dx; 
    p.dz = dz; 
    p.tol = 1e-5;
    p.Tht_amb = Tht_amb;

    p.outfreq = 5;
    p.outvars = {
      {u,   {.name = "u",   .unit = "m/s"}}, 
      {w,   {.name = "w",   .unit = "m/s"}}, 
      {tht, {.name = "tht", .unit = "K"  }}
    };
    p.gnuplot_output = "figure_mr_%s_%d.svg";

    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }

  { // conjugate residual
    using solver_t = output::gnuplot<
      bombel<
	solvers::pressure_cr<
	  solvers::mpdata_2d<real_t, n_iters, n_eqs>, 
	  u, w
        >
      >
    >;
    solver_t::params_t p;

    p.dt = dt; 
    p.dx = dx; 
    p.dz = dz; 
    p.Tht_amb = Tht_amb;
    p.tol = 1e-5;

    p.outfreq = 5;
    p.outvars = {
      {u,   {.name = "u",   .unit = "m/s"}}, 
      {w,   {.name = "w",   .unit = "m/s"}}, 
      {tht, {.name = "tht", .unit = "K"  }}
    };
    p.gnuplot_output = "figure_cr_%s_%d.svg";

    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }

  { // conjugate residual + preconditioner
    using solver_t = output::gnuplot<
      bombel<
        solvers::pressure_pc<
	  solvers::mpdata_2d<real_t, n_iters, n_eqs>,
	  u, w
        >
      >
    >;
    solver_t::params_t p;

    p.dt = dt; 
    p.dx = dx; 
    p.dz = dz; 
    p.Tht_amb = Tht_amb;
    p.tol = 1e-5;
    p.pc_iters = 2;

    p.outfreq = 5;
    p.outvars = {
      {u,   {.name = "u",   .unit = "m/s"}}, 
      {w,   {.name = "w",   .unit = "m/s"}}, 
      {tht, {.name = "tht", .unit = "K"  }}
    };
    p.gnuplot_output = "figure_pc_%s_%d.svg";

    slvs.push_back(new concurr::threads<
      solver_t, 
      bcond::cyclic, // X
      bcond::cyclic  // Y
    >(nx, ny, p));
  }

  for (auto &slv : slvs)
  {
    // initial condition
    {
      blitz::firstIndex i;
      blitz::secondIndex j;

      slv.state(tht) = Tht_amb 
	+ exp( -sqr(i-nx/2.) / (2.*pow(nx/20, 2))
	       -sqr(j-ny/4.) / (2.*pow(ny/20, 2)) )

   	- exp( -sqr(i-nx/3.) / (2.*pow(nx/20, 2))
	       -sqr(j-3*ny/4.) / (2.*pow(ny/20, 2)) )
   	- exp( -sqr(i-2*nx/3.) / (2.*pow(nx/20, 2))
	       -sqr(j-3*ny/4.) / (2.*pow(ny/20, 2)) )

   	+ exp( -sqr(i-2*nx/3.) / (2.*pow(nx/20, 2))
	       -sqr(j-2*ny/4.) / (2.*pow(ny/20, 2)) )
   	+ exp( -sqr(i-2*nx/3.) / (2.*pow(nx/20, 2))
	       -sqr(j-2*ny/4.) / (2.*pow(ny/20, 2)) )

   	+ exp( -sqr(i-nx/2.) / (2.*pow(nx/20, 2))
	       -sqr(j-ny/2.) / (2.*pow(ny/20, 2)) )

   	+ exp( -sqr(i-nx/4.) / (2.*pow(nx/20, 2))
	       -sqr(j-ny/4.) / (2.*pow(ny/20, 2)) )
   	+ exp( -sqr(i-nx/4.) / (2.*pow(nx/20, 2))
	       -sqr(j-3*ny/4.) / (2.*pow(ny/20, 2)) )

   	- exp( -sqr(i-nx/3.) / (2.*pow(nx/20, 2))
	       -sqr(j-ny/2.) / (2.*pow(ny/20, 2)) )
   	- exp( -sqr(i-2*nx/3.) / (2.*pow(nx/20, 2))
	       -sqr(j-ny/2.) / (2.*pow(ny/20, 2)) )
      ;
      slv.state(u) = real_t(0); 
      slv.state(w) = real_t(0); 
    }

    // integration
    slv.advance(nt); 
  }
};
