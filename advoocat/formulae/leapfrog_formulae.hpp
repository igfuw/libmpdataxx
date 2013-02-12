/** @file
 * @copyright University of Warsaw
 * @section LICENSE
 * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */

#pragma once

#include "../idxperm.hpp"

namespace advoocat
{
  namespace formulae
  {
    namespace leapfrog
    {
      using namespace arakawa_c;
      using idxperm::pi;

      /// Implements a leapfrog-type scheme on an Arakawa-C grid.
      /// \f$
      /// \psi^{n+1}_i = \psi^{n-1}_i - C^{n}_{i} \cdot (\psi^{n}_{i+1} - \psi^{n}_{i-1})
      /// \f$
      /// where C is the average Courant number for Arakawa C grid:
      /// \f$
      /// C^{n}_i=0.5\cdot(C^{n}_{i+1/2} + C^{n}_{i-1/2})
      /// \f$
      template <class arr_1d_t>
      void op_1d(
	const arrvec_t<arr_1d_t> &psi, 
	const int n,
	const arr_1d_t &C, 
	const rng_t &i
      ) { 
	psi[n+1](i) = psi[n-1](i) - (C(i+h) + C(i-h)) / 2 * (psi[n](i+1) - psi[n](i-1)); 
      }
    }; // namespace leapfrog
  }; // namespace formulae
}; // namespace advoocat
