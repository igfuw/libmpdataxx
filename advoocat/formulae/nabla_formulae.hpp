/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include "../blitz.hpp"
#include "../idxperm.hpp"
#include "../arakawa_c.hpp"

namespace advoocat
{
  namespace formulae
  {
    namespace nabla_op 
    {
      using idxperm::pi;

      template <int d, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const real_t dx
      ) return_macro(,
	(
	  x(pi<d>(i+1, j)) - 
	  x(pi<d>(i-1, j))
	) / dx / 2.
      )
      
      template <class arg_t, typename real_t>
      inline auto div(
	const arg_t &x,   // x-component of the vector field
	const arg_t &y,   // y-component of the vector field
	const rng_t &i,
	const rng_t &j,
	const real_t hx,
	const real_t hy
      ) return_macro(;,
	(x(i+1,j  ) - x(i-1,j  )) / hx / 2.
	+
	(y(i  ,j+1) - y(i  ,j-1)) / hy / 2.
      )
    }; // namespace nabla_op
  }; // namespace formulae
}; // namespace advoocat
