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
 * @section DERIVATION
 *
 * A system of two 1-dimensional advection equations representing 
 * coupled harmonic oscillators is considered:
 *
 * \f$ \partial_t \psi + \nabla (\vec{u} \psi) =  \omega \phi \f$
 * 
 * \f$ \partial_t \phi + \nabla (\vec{u} \phi) = -\omega \psi \f$
 *
 * Discretisation in time yields:
 *
 * \f$ \frac{\psi^{n+1} - \psi^n}{\Delta t} + A(\psi^n) = \omega \phi^{n+1} \f$
 *
 * \f$ \frac{\phi^{n+1} - \phi^n}{\Delta t} + A(\phi^n) = - \omega \psi^{n+1} \f$
 * 
 * and after some regrouping:
 *
 * \f$ \psi^{n+1} = \Delta t \cdot \omega \phi^{n+1} + \left.\psi^{n+1}\right|_{RHS=0}\f$
 *
 * \f$ \phi^{n+1} = - \Delta t \cdot \omega \psi^{n+1} + \left.\phi^{n+1}\right|_{RHS=0}\f$
 * 
 * solving for \f$ \psi^{n+1} \f$ and \f$ \phi^{n+1} \f$ yields:
 *
 * \f$ \psi^{n+1} = \Delta t \cdot \omega \left( \left.\phi^{n+1}\right|_{RHS=0} - \Delta t \cdot \omega \psi^{n+1} \right) + \left.\psi^{n+1}\right|_{RHS=0} \f$
 *
 * \f$ \phi^{n+1} = - \Delta t \cdot \omega \left( \left.\psi^{n+1}\right|_{RHS=0} + \Delta t \cdot \omega \phi^{n+1} \right) + \left.\phi^{n+1}\right|_{RHS=0}\f$
 *
 * what can be further rearranged to yield:
 *
 * \f$ \psi^{n+1} = \left[ \Delta t \cdot \omega \left.\phi^{n+1}\right|_{RHS=0} + \left.\psi^{n+1}\right|_{RHS=0} \right] / \left[ 1 + \Delta t^2 \cdot \omega^2 \right] \f$
 * 
 * \f$ \phi^{n+1} = \left[ - \Delta t \cdot \omega \left.\psi^{n+1}\right|_{RHS=0} + \left.\phi^{n+1}\right|_{RHS=0} \right] / \left[ 1 + \Delta t^2 \cdot \omega^2 \right] \f$
 *
 * what is represented by the forcings() method in the example below.
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

// auto-deallocating containers
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/assign/ptr_map_inserter.hpp>

// 
#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

enum {psi, phi};
using real_t = double;
using namespace advoocat;

template <class inhomo_solver_t>
class coupled_harmosc : public inhomo_solver_t
{
  using parent_t = inhomo_solver_t;
  using real_t = typename inhomo_solver_t::real_t;
  using arr_1d_t = typename inhomo_solver_t::mem_t::arr_t;

  real_t omega;
  arr_1d_t tmp;

  void forcings(real_t dt)
  {
    auto Psi = this->psi(psi);
    auto Phi = this->psi(phi);

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

  struct params_t : parent_t::params_t { real_t omega; };

  // ctor
  coupled_harmosc(typename parent_t::mem_t &mem, const rng_t &i, params_t p) :
    parent_t(mem, i, p),
    omega(p.omega), 
    tmp(mem.tmp[std::string(__FILE__)][0][0]) 
  {
  }

  static void alloc(
    typename parent_t::mem_t &mem,
    const int nx
  )
  {
    parent_t::alloc(mem, nx);
    const std::string file(__FILE__);
    mem.tmp[file].push_back(new arrvec_t<arr_1d_t>()); 
    mem.tmp[file].back().push_back(new arr_1d_t( rng_t(0, nx-1) )); 
  }
};

int main() 
{
  const int nx = 1000, nt = 750, n_out=10;
  const real_t C = .5, dt = 1;
  const real_t omega = 2*pi<real_t>() / dt / 400;
 
  const int n_iters = 3, n_eqs = 2;

  using solver_t = coupled_harmosc<
    solvers::inhomo_solver_naive< // TODO: plot for both naive and non-naive solver
      solvers::donorcell_1d<
//      solvers::mpdata_1d<
//        n_iters, 
        bcond::cyclic_1d<real_t>, 
        sharedmem_1d<n_eqs, real_t>
      >
    >
  >;

  solver_t::params_t p;
  p.dt = dt;
  p.omega = omega;
  concurr::openmp<solver_t> slv(nx, p);

  Gnuplot gp;
  gp << "set term svg size 1000,500 dynamic enhanced\n" 
     << "set output 'figure.svg'\n";

  gp << "set grid\n";

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
