/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a minimalistic model of a harmonic oscillator
 * (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
 *
 * \image html "../../tests/harmosc/figure_euler_it=1.svg"
 * \image html "../../tests/harmosc/figure_euler_it=1.svg"
 * \image html "../../tests/harmosc/figure_strang_it=2.svg"
 * \image html "../../tests/harmosc/figure_strang_it=2.svg"
 */

#include <libmpdata++/solvers/adv/mpdata_1d.hpp>
#include <libmpdata++/solvers/adv/donorcell_1d.hpp>
#include <libmpdata++/bcond/cyclic_1d.hpp>
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include "coupled_harmosc.hpp"

using real_t = double;
using namespace libmpdataxx;

enum {psi, phi};

template <class T>
void setopts(T &p, std::string name)
{
  p.dt = 1;
  p.omega = 2*pi<real_t>() / p.dt / 400;
  p.gnuplot_output = "figure_" + name + ".svg";
  p.outfreq = 10;
  p.outvars = {
    {psi, {.name = "psi", .unit = "1"}},
    {phi, {.name = "phi", .unit = "1"}}
  };
  p.gnuplot_command = "plot";
}

int main() 
{
  const int nx = 1000, nt = 750;
  const real_t C = .5;
 
  boost::ptr_vector<concurr::any<real_t, 1>> slvs;

  { // euler / donor-cell
    using solver_t = output::gnuplot<coupled_harmosc<real_t, /* n_iters = */ 1, solvers::euler, psi, phi>>;
    solver_t::params_t p; 
    setopts(p, "euler_it=1");
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }
  { // euler / mpdata
    using solver_t = output::gnuplot<coupled_harmosc<real_t, /* n_iters = */ 2, solvers::euler, psi, phi>>;
    solver_t::params_t p; 
    setopts(p, "euler_it=2");
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }
  { // strang / donor-cell
    using solver_t = output::gnuplot<coupled_harmosc<real_t, /* n_iters = */ 1, solvers::strang, psi, phi>>;
    solver_t::params_t p;
    setopts(p, "strang_it=1");
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }
  { // strang / mpdata
    using solver_t = output::gnuplot<coupled_harmosc<real_t, /* n_iters = */ 2, solvers::strang, psi, phi>>;
    solver_t::params_t p; 
    setopts(p, "strang_it=2");
    slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(nx, p));
  }

  for (auto &slv : slvs)
  {
    // initial condition
    {
      blitz::firstIndex i;
      slv.state(psi) = pow(sin((i+.5) * pi<real_t>() / nx), 300);
      slv.state(phi) = real_t(0);
    }
    slv.courant() = C;

    // integration
    slv.advance(nt);
  }
}
