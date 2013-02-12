/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a bombel 
 */

// advection (<> should be used instead of "" in normal usage) 
#include "advoocat/solvers/mpdata_2d.hpp"
#include "advoocat/solvers/solver_inhomo.hpp"
#include "advoocat/solvers/solver_pressure_mr.hpp"
#include "advoocat/solvers/solver_pressure_cr.hpp"
#include "advoocat/solvers/solver_pressure_pc.hpp"
#include "advoocat/bcond/cyclic_2d.hpp"
#include "advoocat/concurr/openmp.hpp"
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
const int n_iters = 2;

using namespace advoocat;

using parent_t_ = 
//solvers::pressure_mr<
solvers::pressure_cr<
//solvers::pressure_pc<
  solvers::inhomo_solver<
    solvers::mpdata_2d<n_iters, sharedmem_2d<3, real_t>>
  >, u, w, tht
>;

// TODO: separate file
class bombel : public parent_t_
{
  using parent_t = parent_t_; // parent_t is a typedef in parent_t_

  real_t Tht_amb;

  void forcings(real_t dt)  //explicit forcings (to be applied before the eliptic solver)
  {
    auto W   = this->psi(w);
    auto Tht = this->psi(tht);

    W += (dt * si:: seconds) * formulae::g<real_t>() * si::seconds / si::metres * (Tht - Tht_amb) / Tht_amb;
  }

  public:

  struct params_t : parent_t::params_t { real_t Tht_amb; };

  // ctor
  bombel(
    parent_t::mem_t &mem, 
    parent_t::bc_p &bcxl,
    parent_t::bc_p &bcxr,
    parent_t::bc_p &bcyl,
    parent_t::bc_p &bcyr,
    const rng_t &i,
    const rng_t &j,
    const params_t &p
  ) :
    parent_t(mem, bcxl, bcxr, bcyl, bcyr, i, j, p),
    Tht_amb(p.Tht_amb)
  {}
};

int main() 
{
  const int nx = 100, ny = 100, nt = 40, n_out=1;
//  const int nx = 50, ny = 50, nt = 41, n_out=10;

  rng_t i(0, nx-1);
  rng_t j(0, ny-1);
  const real_t halo = 1;  

  typename bombel::params_t p;
  p.dt = .1;
  //ambient state (constant thoughout the domain)
  p.Tht_amb = 300;
  //p.Prs_amb = formulae::diagnose::p(p.Tht_amb);
  concurr::openmp<
    bombel,
    bcond::cyclic,
    bcond::cyclic
  > solver(nx, ny, p);

  // initial condition
  {
    blitz::firstIndex i;
    blitz::secondIndex j;

    solver.state(tht) = p.Tht_amb 
      + exp( -sqr(i-nx/2.) / (2.*pow(nx/20, 2))
             -sqr(j-ny/4.) / (2.*pow(ny/20, 2)) )
    ;
    solver.state(u) = real_t(0); 
    solver.state(w) = real_t(0); 
  }

  //ploting
  Gnuplot gp;
  gp << "reset\n"
     << "set term svg size 2000,750 dynamic\n"
     << "set output 'figure.svg'\n"
     << "set multiplot layout 1,3 columnsfirst\n"
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

  std::string binfmt;
  binfmt = gp.binfmt(solver.state());

  // integration
//  for (int t = 1; t <= 5; ++t)
//  {
    solver.advance(nt); // 1 tymczasowo
//   if (t % n_out == 0 /*&& t != 0*/)  
//   {    
//      gp << "set title 'tht @ t=" << t+1 << "'\n"
      gp << "set title 'tht @ t=" << std::setprecision(3) << nt * p.dt << "'\n"
//         << "set cbrange [298.5:302]\n"
         << "splot '-' binary" << binfmt << "with image notitle\n";
      gp.sendBinary(solver.state(tht).copy());
      gp << "set title 'u @ t=" << std::setprecision(3) << nt * p.dt << "'\n"
//         << "set cbrange [-.03:.03]\n"
         << "splot '-' binary" << binfmt << "with image notitle\n";
      gp.sendBinary(solver.state(u).copy());
      gp << "set title 'w @ t=" <<std::setprecision(3) << nt * p.dt << "'\n"
//         << "set cbrange [-.03:.07]\n"
         << "splot '-' binary" << binfmt << "with image notitle\n";
      gp.sendBinary(solver.state(w).copy());
};
