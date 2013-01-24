/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a bombel 
 */

// advection (<> should be used instead of "" in normal usage) 
#include "advoocat/mpdata_2d.hpp"
#include "advoocat/solver_pressure.hpp"
#include "advoocat/cyclic_2d.hpp"
//gradient
#include "advoocat/nabla_formulae.hpp"
//physical constants
#include "advoocat/phc.hpp"
//theta->pressure
#include "advoocat/diagnose_formulae.hpp"

// plotting
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

// auto-deallocating containers
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

enum {u, w, tht};  //eqations
enum {x, z};       //dimensions

template <int n_iters, typename real_t>
using parent = 
pressure_solver<
  inhomo_solver<
    solvers::mpdata_2d<
      n_iters, 
      cyclic_2d<x, real_t>, 
      cyclic_2d<z, real_t>,
      3, 
      real_t
    >
  >, u, w, tht, x, z
>;

template <int n_iters, typename real_t = float>
class bombel : public parent<n_iters, real_t>
{
  using arr_2d_t = typename parent<n_iters, real_t>::arr_t;

  void forcings(real_t dt)
  {
    auto W   = this->state(w);
    auto U   = this->state(u);
    auto Tht = this->state(tht);

    //TODO units, physical constants
    const real_t Tht_amb = 300;     //[K]

    this->xchng(tht); //for gradient

    // diagnose pressure from theta field
    arr_2d_t Prs(this->nx, this->ny+2);
    Prs = diagnose::p(Tht);

    //reference state for pressure
    blitz::Array<real_t, 1> Prs_amb(this->ny+2);  //[Pa]  TODO units, physical constants
    Prs_amb = 81138.2;
    blitz::secondIndex k;

    rng_t l(1, this->ny-1);

    //TODO  can't work yet - no pressure gradient force
    W += (dt * si:: seconds) * phc::g<real_t>() * si::seconds / si::metres * (Tht - Tht_amb) / Tht_amb
        - dt * nabla_op::grad<1>(/*diagnose::p(Tht, this->i, this->j)*/ (Prs - Prs_amb(k) / Prs_amb(k)), l, this->i, real_t(1))
    ;

    //TODO forcings for theta
    // Tht -= dt * (U*nabla_op::grad<0>(Tht_amb) + W*nabla_op::grad<1>(Tht_amb))
  }

  public:

  bombel(int nx, int ny, real_t dt) :
    parent<n_iters, real_t>(nx, ny, dt)
  {
  }
};

int main() 
{
  const int nx = 50, ny = 50, nt = 4, n_out=1;
//  const int nx = 50, ny = 50, nt = 41, n_out=10;
  using real_t = float;
  const real_t dt = .1;

  rng_t i(0, nx-1);
  rng_t j(0, ny-1);
  const real_t halo = 1;  
 
  bombel<2> solver(nx, ny, dt);

  // initial condition
  {
    blitz::firstIndex i;
    blitz::secondIndex j;

    solver.state(tht) = real_t(300) 
      + exp( -sqr(i-nx/2.) / (2.*pow(nx/10, 2))
               -sqr(j-ny/3.) / (2.*pow(ny/10, 2)) )
    ;
    solver.state(u) = real_t(0); 
    solver.state(w) = real_t(0); 
  }

  //ploting
  Gnuplot gp;
  gp << "reset\n"
     << "set term svg size 2000,1000 dynamic\n"
     << "set output 'figure.svg'\n"
     << "set multiplot layout 2,5 columnsfirst\n"
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
  for (int t = 0; t <= nt; ++t)
  {
    solver.solve(1); // 1 tymczasowo

   if (t % n_out == 0 /*&& t != 0*/)  
   {    
      gp << "set title 'tht @ t=" << t+1 << "'\n"
//         << "set cbrange [286:290]\n"
         << "splot '-' binary" << binfmt << "with image notitle\n";
      gp.sendBinary(solver.state(tht).copy());
      gp << "set title 'w @ t=" << t+1 << "'\n"
//         << "set cbrange [-.2:.4]\n"
         << "splot '-' binary" << binfmt << "with image notitle\n";
      gp.sendBinary(solver.state(w).copy());
   }
  }
};
