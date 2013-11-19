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

const int nt = 750;
const real_t C = .5;

template <class solver_t, class slvs_t>
void add_solver(slvs_t &slvs, std::string name, const int n_iters)
{
  typename solver_t::params_t p; 

  // pre instantiation
  p.span = {1000};
  p.n_iters = n_iters;
  p.dt = 1;
  p.omega = 2*pi<real_t>() / p.dt / 400;
  p.gnuplot_output = "figure_" + name + ".svg";
  p.outfreq = 10;
  p.outvars = {
    {psi, {.name = "psi", .unit = "1"}},
    {phi, {.name = "phi", .unit = "1"}}
  };
  p.gnuplot_command = "plot";

  // instantiation
  slvs.push_back(new concurr::threads<solver_t, bcond::cyclic>(p));

  // post instantiation
  {
    blitz::firstIndex i;
    slvs.back().advectee(psi) = pow(sin((i+.5) * pi<real_t>() / p.span[0]), 300);
    slvs.back().advectee(phi) = real_t(0);
  }
  slvs.back().advector() = C;
}

int main() 
{
  boost::ptr_vector<concurr::any<real_t, 1>> slvs;

  add_solver<output::gnuplot<coupled_harmosc<real_t, solvers::euler,  psi, phi>>>(slvs, "euler_it=1",  1);
  add_solver<output::gnuplot<coupled_harmosc<real_t, solvers::euler,  psi, phi>>>(slvs, "euler_it=2",  2);
  add_solver<output::gnuplot<coupled_harmosc<real_t, solvers::strang, psi, phi>>>(slvs, "strang_it=1", 1);
  add_solver<output::gnuplot<coupled_harmosc<real_t, solvers::strang, psi, phi>>>(slvs, "strang_it=2", 2);

  for (auto &slv : slvs) slv.advance(nt);
}
