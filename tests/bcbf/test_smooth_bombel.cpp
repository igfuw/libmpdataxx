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
 */

#include <libmpdata++/solvers/adv/mpdata_2d.hpp>
#include <libmpdata++/solvers/adv+rhs+vip+prs/solver_pressure_mr.hpp>
#include <libmpdata++/solvers/adv+rhs+vip+prs/solver_pressure_cr.hpp>
#include <libmpdata++/solvers/adv+rhs+vip+prs/solver_pressure_pc.hpp>
#include <libmpdata++/bcond/cyclic_2d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

//theta->pressure
#include <libmpdata++/formulae/diagnose_formulae.hpp>

// auto-deallocating containers
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

enum {u, w, tht};  //eqations
enum {x, z};       //dimensions

using real_t = double;
const int n_iters = 2, n_eqs = 3; // TODO: n_eqs should be in bombel!!!

using namespace libmpdataxx;

#include "bombel.hpp"

template <class T>
void setopts(T &p, real_t Tht_amb, std::string name)
{
  p.dt = .1;
  p.dx = p.dz = 1.;
  p.Tht_amb = Tht_amb; 

  p.outfreq = 1;
  p.outvars = {
    {u,   {.name = "u",   .unit = "m/s"}}, 
    {w,   {.name = "w",   .unit = "m/s"}}, 
    {tht, {.name = "tht", .unit = "K"  }}
  };
  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_" + name + "_%s_%d.png";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
  p.gnuplot_contour = true;
  p.gnuplot_term = "png";
// p.gnuplot_cbrange = "[298.5:302]";
}

int main() 
{
  const int nx = 100, ny = 100, nt = 1000;
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
    setopts(p, Tht_amb, "mr");
    p.tol = 1e-5;

    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }
/*
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
    setopts(p, Tht_amb, "cr");
    p.tol = 1e-5;

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
    setopts(p, Tht_amb, "pc");
    p.tol = 1e-5;
    p.pc_iters = 2;

    slvs.push_back(new concurr::threads<
      solver_t, 
      bcond::cyclic, // X
      bcond::cyclic  // Y
    >(nx, ny, p));
  }
*/
  for (auto &slv : slvs)
  {
    // initial condition
    {
      blitz::firstIndex i;
      blitz::secondIndex j;

      slv.state(tht) = Tht_amb 
	+ exp( -sqr(.5+i-nx/2.) / (2.*pow(nx/20, 2))   // TODO: assumed dx=dy=1?
	       -sqr(.5+j-ny/4.) / (2.*pow(ny/20, 2)) )
      ;
      slv.state(u) = real_t(0); 
      slv.state(w) = real_t(0); 
    }

    // integration
    slv.advance(nt); 
  }
};
