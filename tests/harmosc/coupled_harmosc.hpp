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
// TODO: update docs above!

#include <libmpdata++/solvers/mpdata_rhs.hpp>

template <class ct_params_t, int psi, int phi>
class coupled_harmosc : public libmpdataxx::solvers::mpdata_rhs<ct_params_t>
{
  static_assert(ct_params_t::n_eqs == 2, "psi & phi");
  using parent_t = libmpdataxx::solvers::mpdata_rhs<ct_params_t>;
  enum { n = 0 }; // just to make n, n+1 look nice :)

  // member fields
  typename parent_t::real_t omega;
  typename parent_t::arr_t tmp;

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs, 
    typename parent_t::real_t dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at);

    auto Psi = this->state(psi);
    auto Phi = this->state(phi);
    auto &i = this->i;

    switch (at)
    {
      /// explicit solution for R^{n}
      case (n):
      rhs.at(psi)(i) += omega * Phi(i);
      rhs.at(phi)(i) -= omega * Psi(i);
      break;
   
      case (n+1):
      /// implicit solution for R^{n+1}
      /// (consult eq. 28 in Smolarkiewicz 2006, IJNMF)
      rhs.at(psi)(i) += (
	(Psi(i) + dt * omega * Phi(i)) / (1 + pow(dt * omega, 2))
	- 
	Psi(i)
      ) / dt;

      rhs.at(phi)(i) += (
        (Phi(i) - dt * omega * Psi(i)) / (1 + pow(dt * omega, 2))
        -
        Phi(i)
      ) / dt;
      break;
  
      default: assert(false); 
    }
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    typename parent_t::real_t omega; 
  };

  // ctor
  coupled_harmosc(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) :
    parent_t(args, p),
    omega(p.omega)
  {}
};
