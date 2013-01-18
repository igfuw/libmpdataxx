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
#include "advoocat/cyclic.hpp"

// plotting
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

// auto-deallocating containers
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

//TODO as a template
enum {tht, prs, u, w};
enum {x, y};

template <int n_iters>
class bombel : public pressure_solver<solvers::mpdata_2d<n_iters, cyclic_2d<x>, cyclic_2d<y>,4>>
{
  void forcings(real_t dt)
  {
  auto W = this->state(w);
  auto Tht = this->state(tht);
  
  //TODO units, physical constants
  const real_t g = 9.81;         //[m/s]
  const real_t Tht_amb = 287;    //[K]

  //TODO  can't work yet - no pressure gradient force
  W += dt * g * (Tht - Tht_amb) / Tht_amb; 

  //TODO forcings for theta
 
  }

  public:

  bombel(int nx, int ny, real_t dt) :
    pressure_solver<solvers::mpdata_2d<n_iters, cyclic_2d<x>, cyclic_2d<y>, 4>>(nx, ny, dt)
  {
  }
};

int main() 
{
  const int nx = 10, ny = 10, nt = 2, n_out=1;
  const real_t C = .5, dt = 1;

  rng_t i(0, nx-1);
  rng_t j(0, ny-1);
  const real_t halo = 1;  
 
  bombel<2> solver(nx, ny, dt);

  // initial condition
  {
    blitz::firstIndex i;
    blitz::secondIndex j;

    solver.state(tht) = real_t(287) + 2*exp(
      -sqr(i-nx/2.) / (2.*pow(nx/10, 2))
      -sqr(j-ny/4.) / (2.*pow(ny/10, 2))
    );
    solver.state(prs) = real_t(101300);
    solver.state(u) = real_t(0);
    solver.state(w) = real_t(0);
  }

  //ploting
  Gnuplot gp;
  gp << "reset\n"
     << "set term svg size 1000,1000 dynamic\n"
     << "set output 'figure.svg'\n"
     << "set multiplot layout 2,2\n"
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

  // plot initial conditions
  gp << "set title 'theta'\n"
     << "set cbrange [287:289]\n"
     << "splot '-' binary" << binfmt << "with image notitle\n";
  gp.sendBinary(solver.state(tht).copy());

  gp << "set title 'pressure'\n"
     << "set cbrange [101200:101400]\n"
     << "splot '-' binary" << binfmt << "with image notitle\n";
  gp.sendBinary(solver.state(prs).copy());

  // integration
  for (int t = 0; t <= nt; ++t)
  {
    // TODO: trzeba pamiętać o odpowiedniku fill_halos dla courantów

    // uwaga: aktualnie courant() (w przeciwienstwie do state()) 
    // zwraca cala tablice razem z halo!
    solver.courant(x) = real_t(0); 
    solver.courant(y) = real_t(0);

    solver.solve(1); // 1 tymczasowo

    if (t % n_out == 0 && t != 0)  
    {    
      gp << "set cbrange [287:289]\n"
         << "set title ' '\n"
         << "splot '-' binary" << binfmt << "with image notitle\n";
      gp.sendBinary(solver.state(tht).copy());

      gp << "set cbrange [101200:101400]\n"
         << "set title ' '\n"
         << "splot '-' binary" << binfmt << "with image notitle\n";
      gp.sendBinary(solver.state(prs).copy());
    }
  }
};
