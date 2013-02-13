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
 */

#include "advoocat/blitz.hpp"

using namespace advoocat;

template <class inhomo_solver_t, int psi, int phi>
class coupled_harmosc : public inhomo_solver_t
{
  using parent_t = inhomo_solver_t;

  typename parent_t::real_t omega;
  typename parent_t::arr_t tmp;

  void forcings(typename parent_t::real_t dt)
  {
    auto Psi = parent_t::psi(psi);
    auto Phi = parent_t::psi(phi);

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

  struct params_t : parent_t::params_t { typename parent_t::real_t omega; };

  // ctor
  coupled_harmosc(
    typename parent_t::mem_t *mem, 
    typename parent_t::bc_p &bcxl,
    typename parent_t::bc_p &bcxr,
    const rng_t &i, 
    params_t p
  ) :
    parent_t(mem, bcxl, bcxr, i, p),
    omega(p.omega), 
    tmp(mem->tmp[std::string(__FILE__)][0][0]) 
  {}

  static void alloc(
    typename parent_t::mem_t *mem,
    const int nx
  )
  {
    parent_t::alloc(mem, nx);
    const std::string file(__FILE__);
    mem->tmp[file].push_back(new arrvec_t<typename parent_t::arr_t>()); 
    mem->tmp[file].back().push_back(new typename parent_t::arr_t( rng_t(0, nx-1) )); 
  }
};
