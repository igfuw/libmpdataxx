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

#include <libmpdata++/solvers/adv/donorcell_2d.hpp>
#include <libmpdata++/solvers/adv/mpdata_2d.hpp>
#include <libmpdata++/solvers/adv/mpdata_fct_2d.hpp>
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
const int n_iters = 2;

using namespace libmpdataxx;

#include "bombel.hpp"

template <class T>
void setopts(T &p, real_t Tht_amb, std::string name)
{
  p.dt = .5;
  p.dx = p.dz = 10.;
  p.Tht_amb = Tht_amb; 

  p.outfreq = 1;
  p.outvars = {
//    {u,   {.name = "u",   .unit = "m/s"}}, 
//    {w,   {.name = "w",   .unit = "m/s"}}, 
    {tht, {.name = "tht", .unit = "K"  }}
  };
  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_" + name + "_%s_%d.png";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
  p.gnuplot_contour = true;
//  p.gnuplot_cbrange = "[299:301.5]";
  p.gnuplot_term = "png";
}

int main() 
{
  const int nx = 200, ny = 200, nt = 1212;
  real_t Tht_amb = 300; // ambient state (constant thoughout the domain)

  boost::ptr_vector<concurr::any<real_t, 2>> slvs;

  { // minimum residual
    using solver_t = output::gnuplot<
      bombel<
	solvers::pressure_cr<
//        solvers::donorcell_2d<real_t>,
	  solvers::mpdata_fct_2d<real_t, n_iters>,
	  u, w
        >
      >
    >;
    solver_t::params_t p; 
    setopts(p, Tht_amb, "cr");
    p.tol = 1e-5;

    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }
  for (auto &slv : slvs)
  {
    // initial condition
    {
      blitz::firstIndex i;
      blitz::secondIndex j;

      slv.advectee(tht) = Tht_amb 
	+ .5 * exp( -sqr(.5+i-nx/2.) / (2.*pow(nx/20, 2))   // TODO: assumed dx=dy=1?
	       	    -sqr(.5+j-ny/4.) / (2.*pow(ny/20, 2)) )
      ;
      slv.advectee(u) = real_t(0); 
      slv.advectee(w) = real_t(0); 
    }

    // integration
    slv.advance(nt); 
  }
};
