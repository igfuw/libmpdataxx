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
#include "advoocat/solvers/solver_inhomo.hpp"
#include "advoocat/solvers/solver_pressure_mr.hpp"
#include "advoocat/solvers/solver_pressure_cr.hpp"
#include "advoocat/solvers/solver_pressure_pc.hpp"
#include "advoocat/bcond/cyclic_2d.hpp"
#include "advoocat/concurr/threads.hpp"
//gradient
#include "advoocat/formulae/nabla_formulae.hpp"
//theta->pressure
#include "advoocat/formulae/diagnose_formulae.hpp"
//physical constants
#include "advoocat/formulae/phc.hpp"

// plotting
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

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
  const int nx = 100, ny = 100, nt = 40, n_out=1;
//  const int nx = 20, ny = 20, nt = 1, n_out=1;

  rng_t i(0, nx-1);
  rng_t j(0, ny-1);
  const real_t halo = 1;  

  // TODO p.Prs_amb = formulae::diagnose::p(p.Tht_amb);
  real_t dt = .1; // timestep
  real_t Tht_amb = 300; // ambient state (constant thoughout the domain)

  boost::ptr_vector<concurr::any<real_t, 2>> slvs;

  { // minimum residual
    using solver_t = bombel<
      solvers::pressure_mr<
	solvers::inhomo_solver<
	  solvers::mpdata_2d<real_t, n_iters, n_eqs>, solvers::strang
	>, u, w, tht
      >
    >;
    solver_t::params_t p; p.dt = dt; p.Tht_amb = Tht_amb;
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }
  { // conjugate residual
    using solver_t = bombel<
      solvers::pressure_cr<
	solvers::inhomo_solver<
	  solvers::mpdata_2d<real_t, n_iters, n_eqs>, solvers::strang
	>, u, w, tht
      >
    >;
    solver_t::params_t p; p.dt = dt; p.Tht_amb = Tht_amb;
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }

  { // conjugate residual + preconditioner
    using solver_t = bombel<
      solvers::pressure_pc<
	solvers::inhomo_solver<
	  solvers::mpdata_2d<real_t, n_iters, n_eqs>, solvers::strang
	>, u, w, tht
      >
    >;
    solver_t::params_t p; p.dt = dt; p.Tht_amb = Tht_amb;
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic, bcond::cyclic>(nx, ny, p));
  }

  //ploting
  Gnuplot gp;
  gp << "reset\n"
     << "set term svg size 2000,750 dynamic\n"
     << "set output 'figure.svg'\n"
     << "set multiplot layout " << slvs.size() << ",3\n" // columnsfirst\n"
     << "set grid\n"
     << "set xlabel 'X'\n"
     << "set ylabel 'Y'\n"
     << "set xrange [0:" << nx-1 << "]\n"
     << "set yrange [0:" << ny-1 << "]\n"
     // progressive-rock connoisseur palette ;)
     << "set palette defined (0 '#ffffff', 1 '#993399', 2 '#00CCFF', 3 '#66CC00', 4 '#FFFF00', 5 '#FC8727', 6 '#FD0000')\n" 
     << "set view map\n"
     << "set key font \",5\"\n "
     << "set contour base\n" 
     << "set nosurface\n"
     << "set cntrparam levels 0\n";

  for (auto &slv : slvs)
  {
    // initial condition
    {
      blitz::firstIndex i;
      blitz::secondIndex j;

      slv.state(tht) = Tht_amb 
	+ exp( -sqr(i-nx/2.) / (2.*pow(nx/20, 2))
	       -sqr(j-ny/4.) / (2.*pow(ny/20, 2)) )
      ;
      slv.state(u) = real_t(0); 
      slv.state(w) = real_t(0); 
    }

    std::string binfmt;
    binfmt = gp.binfmt(slv.state());

    // integration
  //  for (int t = 1; t <= 5; ++t)
  //  {
      slv.advance(nt); // 1 tymczasowo
  //   if (t % n_out == 0 /*&& t != 0*/)  
  //   {    
  //      gp << "set title 'tht @ t=" << t+1 << "'\n"
	gp << "set title 'tht @ t=" << std::setprecision(3) << nt * dt << "'\n"
  //         << "set cbrange [298.5:302]\n"
	   << "splot '-' binary" << binfmt << "with image notitle\n";
	gp.sendBinary(slv.state(tht).copy());
	gp << "set title 'u @ t=" << std::setprecision(3) << nt * dt << "'\n"
  //         << "set cbrange [-.03:.03]\n"
	   << "splot '-' binary" << binfmt << "with image notitle\n";
	gp.sendBinary(slv.state(u).copy());
	gp << "set title 'w @ t=" <<std::setprecision(3) << nt * dt << "'\n"
  //         << "set cbrange [-.03:.07]\n"
	   << "splot '-' binary" << binfmt << "with image notitle\n";
	gp.sendBinary(slv.state(w).copy());
  }
};
