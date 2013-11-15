/** 
 * @file
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
 */

#include <libmpdata++/solvers/adv+rhs/solver_inhomo.hpp>
#include <libmpdata++/blitz.hpp>

using namespace libmpdataxx;

template <typename real_t, int n_iters, solvers::inhomo_e inhomo, int psi, int phi>
class coupled_harmosc : public solvers::inhomo_solver<solvers::mpdata_1d<real_t, n_iters>, inhomo>
{
  using parent_t = solvers::inhomo_solver<solvers::mpdata_1d<real_t, n_iters>, inhomo>;

  typename parent_t::real_t omega;
  typename parent_t::arr_t tmp;

  void apply_forcings(typename parent_t::real_t dt)
  {
    auto Psi = this->state(psi);
    auto Phi = this->state(phi);
    rng_t &i = this->i;

    tmp(i) = Psi(i);
    ///   (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
    // explicit part
    Psi(i) += dt * omega * Phi(i);
    // implicit part
    Psi(i) /= (1 + pow(dt * omega, 2));

    // explicit part
    Phi(i) += - dt * omega * tmp(i);
    // implicit part
    Phi(i) /= (1 + pow(dt * omega, 2));
  }

  public:

  struct params_t : parent_t::params_t 
  { 
    typename parent_t::real_t omega; 

    // ctor setting n_eqs
    params_t() { this->n_eqs = 2; }
  };

  // ctor
  coupled_harmosc(
    typename parent_t::ctor_args_t args,
    const params_t &p
  ) :
    parent_t(args, p),
    omega(p.omega), 
    tmp(args.mem->tmp[__FILE__][0][0]) 
  {
    assert(p.n_eqs == 2);
  }

  static void alloc(
    typename parent_t::mem_t *mem,
    const int nx
  )
  {
    // TODO: move to inhomo!
    parent_t::alloc(mem, nx);
    mem->tmp[__FILE__].push_back(new arrvec_t<typename parent_t::arr_t>()); 
    mem->tmp[__FILE__].back().push_back(new typename parent_t::arr_t( rng_t(0, nx-1) )); 
  }
};
