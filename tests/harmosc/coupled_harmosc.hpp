/** 
 * @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *
 * @brief a minimalistic model of a harmonic oscillator
 * (consult eq. 28-30 in Smolarkiewicz 2006, IJNMF)
 *
 */
//

//<listing-1>
#include <libmpdata++/solvers/mpdata_rhs.hpp>

template <class ct_params_t>
struct coupled_harmosc : public 
  libmpdataxx::solvers::mpdata_rhs<ct_params_t>
{
  // aliases
  using parent_t = 
    libmpdataxx::solvers::mpdata_rhs<ct_params_t>;
  using ix = typename ct_params_t::ix;

  // member fields
  typename ct_params_t::real_t omega;

  // method called by mpdata_rhs
  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs, 
    const typename parent_t::real_t &dt,
    const int &at
  ) {
    parent_t::update_rhs(rhs, dt, at);

    // just to shorten code
    auto psi = this->psi_n(ix::psi);
    auto phi = this->psi_n(ix::phi);
    auto &i = this->i;

    switch (at) {
      case (0): // explicit solution for R^{n}
      rhs.at(ix::psi)(i) += omega * phi(i);
      rhs.at(ix::phi)(i) -= omega * psi(i);
      break;
   
      case (1): // implicit solution for R^{n+1}
      rhs.at(ix::psi)(i) += (
	(psi(i) + dt * omega * phi(i)) 
        / (1 + pow(dt * omega, 2))
	- psi(i)
      ) / dt;

      rhs.at(ix::phi)(i) += (
        (phi(i) - dt * omega * psi(i)) 
        / (1 + pow(dt * omega, 2))
        - phi(i)
      ) / dt;
      break;
    }
  }

  // run-time parameters
  struct rt_params_t : parent_t::rt_params_t { 
    typename ct_params_t::real_t omega = 0; 
  };

  // ctor
  coupled_harmosc(
    typename parent_t::ctor_args_t args,
    const rt_params_t &p
  ) : parent_t(args, p), omega(p.omega)
  { assert(omega != 0); }
};
//</listing-1>
