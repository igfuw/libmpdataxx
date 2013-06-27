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
const int n_iters = 2, n_eqs = 3; // TODO: n_eqs should be in bombel!!!

using namespace libmpdataxx;

#include "bombel.hpp"

template <class T>
void setopts(T &p, real_t Tht_amb, std::string name)
{
  p.dt = .5;
  p.dx = p.dz = 10.; 
  p.Tht_amb = Tht_amb; 

  p.outfreq = 10;
  p.outvars = {
//    {u,   {.name = "u",   .unit = "m/s"}}, 
//    {w,   {.name = "w",   .unit = "m/s"}}, 
    {tht, {.name = "tht", .unit = "K"  }}
  };
  p.gnuplot_view = "map";
  p.gnuplot_output = "figure_" + name + "_%s_%04d.png";
  p.gnuplot_with = "lines";
  p.gnuplot_surface = false;
  p.gnuplot_contour = true;
  p.gnuplot_cbrange = "[299.8:300.6]";
  p.gnuplot_maxcolors = 8;
  p.gnuplot_cntrparam = "levels incremental 299.8, 0.1, 300.6";
  p.gnuplot_term = "png";
}

int main() 
{
  const int r0 = 250; 
  const int nx = 200, ny = 200, nt = 1212;
  real_t Tht_amb = 300; // ambient state (constant thoughout the domain)

  boost::ptr_vector<libmpdataxx::concurr::any<real_t, 2>> slvs;

  // conjugate residual
  using solver_t = output::gnuplot<
    bombel<
      solvers::pressure_cr<
        solvers::mpdata_fct_2d<real_t, n_iters, n_eqs>, 
        //solvers::donorcell_2d<real_t, n_eqs>, 
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

    slv.state(tht) = Tht_amb + where(pow(i * p.dx - 4 * r0 , 2) + pow(j * p.dz - 1.04 * r0 , 2) <= pow(r0, 2) , .5, 0);
    slv.state(u) = real_t(0); 
    slv.state(w) = real_t(0); 

    // integration
    slv.advance(nt); 
  }
};
