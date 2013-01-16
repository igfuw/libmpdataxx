/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a minimalistic model of a harmonic oscillator
 * (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
 *
 * \include "harmosc/test_harmosc.cpp"
 * \image html "../../tests/harmosc/figure.svg"
 */

// advection (<> should be used instead of "" in normal usage)
#include "advoocat/lib.hpp"
#include "advoocat/solver_inhomo.hpp"

// plotting
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

// auto-deallocating containers
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

// 
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

enum {psi, phi};

template <int n_iters>
class coupled_harmosc : public inhomo_solver<solvers::mpdata_1d<n_iters, cyclic_1d, 2>>
{
  real_t omega;
  arr_1d_t tmp;

  void forcings(real_t dt)
  {
    auto Psi = this->state(psi);
    auto Phi = this->state(phi);

    tmp = Psi;
    ///   (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
    // explicit part
    Psi += dt * omega * Phi;
    // implicit part
    Psi /= (1 + pow(dt * omega, 2));

    // explicit part
    Phi += - dt * omega * tmp;
    // implicit part
    Phi /= (1 + pow(dt * omega, 2));
  }

  public:

  coupled_harmosc(int n, real_t dt, real_t omega) :
    inhomo_solver<solvers::mpdata_1d<n_iters, cyclic_1d, 2>>(n, dt),
    omega(omega), tmp(this->state(0).extent(0))
  {
  }
};

int main() 
{
  const int nx = 1000, nt = 750, n_out=10;
  const real_t C = .5, dt = 1;
  const real_t omega = 2*pi<real_t>() / dt / 400;
 
  coupled_harmosc<3> solver(nx, dt, omega);

  Gnuplot gp;
  gp << "set term svg size 1000,500 dynamic enhanced\n" 
     << "set output 'figure.svg'\n";

  gp << "set grid\n";

  // initial condition
  {
    blitz::firstIndex i;
    solver.state(psi) = pow(sin(i * pi<real_t>() / nx), 300);
    solver.state(phi) = real_t(0);
  }
  solver.courant() = C;

  gp << "plot"
     << "'-' lt 1 with lines title 'psi',"
     << "'-' lt 2 with lines title 'phi',"
     << "'-' lt 3 with lines title 'psi^2 + phi^2 + 1'";
  for (int t = n_out; t <= nt; t+=n_out) 
    gp << ", '-' lt 1 with lines notitle"
       << ", '-' lt 2 with lines notitle"
       << ", '-' lt 3 with lines notitle";
  gp << "\n";

  arr_1d_t en(nx);

  // sending initial condition
  gp.send(solver.state(psi));
  gp.send(solver.state(phi));
  en = 1 + pow(solver.state(psi),2) + pow(solver.state(phi),2);
  gp.send(en);

  // integration
  for (int t = n_out; t <= nt; t+=n_out)
  {
    solver.solve(n_out);
    gp.send(solver.state(psi));
    gp.send(solver.state(phi));
    en = 1 + pow(solver.state(psi),2) + pow(solver.state(phi),2);
    gp.send(en);
  }
}
