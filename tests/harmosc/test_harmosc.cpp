/** 
 * @file
 * @example harmosc/test_harmosc.cpp
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a minimalistic model of a harmonic oscillator
 * (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
 *
 * @section FIGURE
 *
 * \image html "../../tests/harmosc/figure.svg"
 */

// advection (<> should be used instead of "" in normal usage)
#include "advoocat/solvers/mpdata_1d.hpp"
#include "advoocat/solvers/donorcell_1d.hpp"
#include "advoocat/solvers/solver_inhomo.hpp"
#include "advoocat/bcond/cyclic_1d.hpp"
#include "advoocat/concurr/openmp.hpp"

// plotting
#define GNUPLOT_ENABLE_BLITZ
#include <gnuplot-iostream/gnuplot-iostream.h>

// 
#include <boost/math/constants/constants.hpp>

#include "coupled_harmosc.hpp"

int main() 
{
  using real_t = double;
  using namespace advoocat;
  using boost::math::constants::pi;

  const int nx = 1000, nt = 750, n_out=10;
  const real_t C = .5, dt = 1;
  const real_t omega = 2*pi<real_t>() / dt / 400;
 
  const int n_iters = 3, n_eqs = 2;

  using solver_t = coupled_harmosc<
    solvers::inhomo_solver_naive< // TODO: plot for both naive and non-naive solver
      solvers::donorcell_1d<sharedmem_1d<n_eqs, real_t>>
    >
  >;

  solver_t::params_t p;
  p.dt = dt;
  p.omega = omega;
  concurr::openmp<solver_t, bcond::cyclic> slv(nx, p);

  Gnuplot gp;
  gp << "set term svg size 1000,500 dynamic enhanced\n" 
     << "set output 'figure.svg'\n";

  gp << "set grid\n";

  // initial condition
  {
    blitz::firstIndex i;
    slv.state(solver_t::psi) = pow(sin(i * pi<real_t>() / nx), 300);
    slv.state(solver_t::phi) = real_t(0);
  }
  slv.courant() = C;

  gp << "plot"
     << "'-' lt 1 with lines title 'psi',"
     << "'-' lt 2 with lines title 'phi',"
     << "'-' lt 3 with lines title 'psi^2 + phi^2 + 1'";
  for (int t = n_out; t <= nt; t+=n_out) 
    gp << ", '-' lt 1 with lines notitle"
       << ", '-' lt 2 with lines notitle"
       << ", '-' lt 3 with lines notitle";
  gp << "\n";

  decltype(slv.state()) en(nx);

  // sending initial condition
  gp.send(slv.state(solver_t::psi));
  gp.send(slv.state(solver_t::phi));
  en = 1 + pow(slv.state(solver_t::psi),2) + pow(slv.state(solver_t::phi),2);
  gp.send(en);

  // integration
  for (int t = n_out; t <= nt; t+=n_out)
  {
    slv.advance(n_out);
    gp.send(slv.state(solver_t::psi));
    gp.send(slv.state(solver_t::phi));
    en = 1 + pow(slv.state(solver_t::psi),2) + pow(slv.state(solver_t::phi),2);
    gp.send(en);
  }
}
