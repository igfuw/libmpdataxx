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

#include "coupled_harmosc.hpp"
#include <libmpdata++/concurr/threads.hpp>
#include <libmpdata++/output/gnuplot.hpp>
using namespace libmpdataxx;

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

using real_t = double;

enum {psi, phi};

const int nt = 750;
const real_t C = .5;

template <
  formulae::opts::opts_t opt, 
  solvers::rhs_scheme_t scheme, 
  class slvs_t
> void add_solver(slvs_t &slvs, std::string name, const int n_iters)
{
  struct ct_params_t : ct_params_default_t
  {
    using real_t = real_t;
    enum { n_dims = 1 };
    enum { n_eqs = 2 };
    enum { opts = opt };
    enum { rhs_scheme = scheme };
  };

  using solver_t = output::gnuplot<coupled_harmosc<ct_params_t, psi, phi>>;
  typename solver_t::rt_params_t p; 

  // run-time parameters
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

  add_solver<formulae::opts::abs, solvers::rhs_scheme_t::euler_b>(slvs, "euler_it=1",  1);
  add_solver<formulae::opts::abs, solvers::rhs_scheme_t::euler_b>(slvs, "euler_it=2",  2);
  add_solver<formulae::opts::abs, solvers::rhs_scheme_t::strang>(slvs, "strang_it=1", 1);
  add_solver<formulae::opts::abs, solvers::rhs_scheme_t::strang>(slvs, "strang_it=2", 2);

  for (auto &slv : slvs) slv.advance(nt);
}
