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
#include "advoocat/bcond/cyclic_1d.hpp"
#include "advoocat/concurr/threads.hpp"

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
 
  enum {psi, phi};

  boost::ptr_vector<concurr::any<real_t, 1>> slvs;

  { // euler / donor-cell
    using solver_t = coupled_harmosc<real_t, /* n_iters = */ 1, solvers::euler, psi, phi>;
    solver_t::params_t p; p.dt = dt; p.omega = omega;
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }
  { // euler / mpdata
    using solver_t = coupled_harmosc<real_t, /* n_iters = */ 2, solvers::euler, psi, phi>;
    solver_t::params_t p; p.dt = dt; p.omega = omega;
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }
  { // strang / donor-cell
    using solver_t = coupled_harmosc<real_t, /* n_iters = */ 1, solvers::strang, psi, phi>;
    solver_t::params_t p; p.dt = dt; p.omega = omega;
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }
  { // strang / mpdata
    using solver_t = coupled_harmosc<real_t, /* n_iters = */ 2, solvers::strang, psi, phi>;
    solver_t::params_t p; p.dt = dt; p.omega = omega;
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }

  Gnuplot gp;
  gp << "set term svg size 1000,500 dynamic enhanced\n" 
     << "set output 'figure.svg'\n"
     << "set grid\n"
     << "set multiplot layout 2,2\n";

  for (auto &slv : slvs)
  {
    // initial condition
    {
      blitz::firstIndex i;
      slv.state(psi) = pow(sin(i * pi<real_t>() / nx), 300);
      slv.state(phi) = real_t(0);
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
    gp.send(slv.state(psi));
    gp.send(slv.state(phi));
    en = 1 + pow(slv.state(psi),2) + pow(slv.state(phi),2);
    gp.send(en);

    // integration
    for (int t = n_out; t <= nt; t+=n_out)
    {
      slv.advance(n_out);
      gp.send(slv.state(psi));
      gp.send(slv.state(phi));
      en = 1 + pow(slv.state(psi),2) + pow(slv.state(phi),2);
      gp.send(en);
    }
  }
}
