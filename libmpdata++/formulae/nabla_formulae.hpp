/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/blitz.hpp>
#include <libmpdata++/formulae/idxperm.hpp>
#include <libmpdata++/formulae/arakawa_c.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace nabla
    {
      using idxperm::pi;

      template <class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const real_t dx
      ) return_macro(,
	(
	  x(i+1) - 
	  x(i-1)
	) / dx / 2.
      )

      // 2D version
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
      
      // 3D version
      template <int d, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k,
	const real_t dx
      ) return_macro(,
	(
	  x(pi<d>(i+1, j, k)) - 
	  x(pi<d>(i-1, j, k))
	) / dx / 2.
      )
      
      // 2D version
      template <class arg_t, typename real_t>
      inline auto div(
	const arg_t &x,   // x-component of the vector field
	const arg_t &y,   // y-component of the vector field
	const rng_t &i,
	const rng_t &j,
	const real_t hx,
	const real_t hy
      ) return_macro(,
	(x(i+1,j  ) - x(i-1,j  )) / hx / 2.
	+
	(y(i  ,j+1) - y(i  ,j-1)) / hy / 2.
      )
      
      // 3D version
      template <class arg_t, typename real_t>
      inline auto div(
	const arg_t &x,   // x-component of the vector field
	const arg_t &y,   // y-component of the vector field
	const arg_t &z,   // z-component of the vector field
	const rng_t &i,
	const rng_t &j,
	const rng_t &k,
	const real_t hx,
	const real_t hy,
	const real_t hz
      ) return_macro(,
	(x(i+1, j, k) - x(i-1, j, k)) / hx / 2.
	+
	(y(i, j+1, k) - y(i, j-1, k)) / hy / 2.
	+
	(z(i, j, k+1) - z(i, j, k-1)) / hz / 2.
      )
    } // namespace nabla_op
  } // namespace formulae
} // namespace libmpdataxx
