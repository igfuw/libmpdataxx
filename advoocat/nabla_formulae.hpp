/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "blitz.hpp"
#include "pi.hpp"
#include "arakawa_c.hpp"

namespace nabla_op 
{
  template <int d, class arg_t, typename real_t>
  inline auto grad(
    const arg_t &x,
    const rng_t &i,
    const rng_t &j,
    const real_t h
  ) return_macro(
    (
      x(pi<d>(i+1, j)) - 
      x(pi<d>(i-1, j))
    ) / h / 2.
  )
}; //nabla

