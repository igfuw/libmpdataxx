/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief test for pressure solvers (as in Smolarkiewicz & Pudykiewicz 1992, fig.3) 
 * buoyant convection in Boussinesq flow
 */

//#include <libmpdata++/solvers/adv/donorcell_2d.hpp>
//#include <libmpdata++/solvers/adv/mpdata_2d.hpp>
#include <libmpdata++/solvers/adv/mpdata_fct_2d.hpp>
//#include <libmpdata++/solvers/adv+rhs+vip+prs/solver_pressure_mr.hpp>
#include <libmpdata++/solvers/adv+rhs+vip+prs/solver_pressure_cr.hpp>
//#include <libmpdata++/solvers/adv+rhs+vip+prs/solver_pressure_pc.hpp>
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

using namespace libmpdataxx;

#include "bombel.hpp"

template <class T>
void setopts(T &p, real_t Tht_amb, std::string name)
{
  p.n_iters = 2;

  p.dt = .75;
  p.dx = p.dz = 10.; 
  p.Tht_amb = Tht_amb; 

  p.outfreq = 100; //12;
  p.outvars = {
//    {u,   {.name = "u",   .unit = "m/s"}}, 
//    {w,   {.name = "w",   .unit = "m/s"}}, 
    {tht, {.name = "tht", .unit = "K"  }}
  };
  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_" + name + "_%s_%04d.svg";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
  p.gnuplot_contour = true;
  p.gnuplot_cbrange = "[299.85 : 300.65]";
//  p.gnuplot_cbrange = "[299.85 - 299 : 300.65 - 299]";
  p.gnuplot_maxcolors = 8;
  p.gnuplot_cntrparam = "levels incremental 299.85, 0.1, 300.65";
//  p.gnuplot_cntrparam = "levels incremental 299.85 - 299, 0.1, 300.65 - 299";
  p.gnuplot_term = "svg";
}

int main() 
{
  const int r0 = 250; 
  const int nx = 200, ny = 200, nt = 800;
  real_t Tht_amb = 300; //1; //300; // ambient state (constant thoughout the domain)

  boost::ptr_vector<libmpdataxx::concurr::any<real_t, 2>> slvs;

  // conjugate residual
  using solver_t = output::gnuplot<
    bombel<
      solvers::pressure_cr<
        solvers::mpdata_fct_2d<real_t, formulae::opts::iga>, 
        //solvers::donorcell_2d<real_t>, 
	u, w
      >
    >
  >;
  solver_t::params_t p;
  setopts(p, Tht_amb, "cr");
  p.tol = 1e-5;

  slvs.push_back(new libmpdataxx::concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));

  for (auto &slv : slvs)
  {
    // initial condition
    blitz::firstIndex i;
    blitz::secondIndex j;

    slv.advectee(tht) = Tht_amb + where(
      // if
      pow((i+.5) * p.dx - 4 * r0 , 2) + pow((j+.5) * p.dz - 1.04 * r0 , 2) <= pow(r0, 2), 
      // then
      .5, 
      // else
      0
    );
    slv.advectee(u) = real_t(0); 
    slv.advectee(w) = real_t(0); 

    // integration
    slv.advance(nt); 
  }
};
