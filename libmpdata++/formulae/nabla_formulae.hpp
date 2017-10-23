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
      using arakawa_c::h;

      // 1D version
      template <int ord = 2, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const real_t dx,
        typename std::enable_if<ord == 2>::type* = 0
      ) return_macro(,
	(
	  x(i+1) - 
	  x(i-1)
	) / dx / 2.
      )
      
      // 4th order version of the above
      template <int ord, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const real_t dx,
        typename std::enable_if<ord == 4>::type* = 0
      ) return_macro(,
	(
	  -   x(i+2) +
	  8 * x(i+1) - 
	  8 * x(i-1) +
	      x(i-2)
	) / (12 * dx)
      )

      // 2D version
      template <int d, int ord = 2, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const real_t dx,
        typename std::enable_if<ord == 2>::type* = 0
      ) return_macro(,
	(
	  x(pi<d>(i+1, j)) - 
	  x(pi<d>(i-1, j))
	) / dx / 2.
      )

      // 4th order version of the above
      template <int d, int ord, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const real_t dx,
        typename std::enable_if<ord == 4>::type* = 0
      ) return_macro(,
	(
	  -   x(pi<d>(i+2, j)) + 
	  8 * x(pi<d>(i+1, j)) - 
	  8 * x(pi<d>(i-1, j)) +
	      x(pi<d>(i-2, j))
	) / (12 * dx)
      )
      
      // 3D version
      template <int d, int ord = 2, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k,
	const real_t dx,
        typename std::enable_if<ord == 2>::type* = 0
      ) return_macro(,
	(
	  x(pi<d>(i+1, j, k)) - 
	  x(pi<d>(i-1, j, k))
	) / dx / 2.
      )
      
      // 4th order version of the above
      template <int d, int ord, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k,
	const real_t dx,
        typename std::enable_if<ord == 4>::type* = 0
      ) return_macro(,
	(
	  -   x(pi<d>(i+2, j, k)) +
	  8 * x(pi<d>(i+1, j, k)) - 
	  8 * x(pi<d>(i-1, j, k)) +
	      x(pi<d>(i-2, j, k))
	) / (12 * dx)
      )
      
      // compact gradients

      template <class arg_t, typename real_t>
      inline auto grad_cmpct(
	const arg_t &x,
	const rng_t &i,
	const real_t dx
      ) return_macro(,
	(
	  x(i+1) - 
	  x(i)
	) / dx
      )

      // 2D version
      template <int d, class arg_t, typename real_t>
      inline auto grad_cmpct(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const real_t dx
      ) return_macro(,
	(
	  x(pi<d>(i+1, j)) - 
	  x(pi<d>(i  , j))
	) / dx
      )
      
      // 3D version
      template <int d, class arg_t, typename real_t>
      inline auto grad_cmpct(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k,
	const real_t dx
      ) return_macro(,
	(
	  x(pi<d>(i+1, j, k)) - 
	  x(pi<d>(i, j, k))
	) / dx
      )
      
      // helper function to calculate gradient components of a scalar field
      
      // 1D version
      template <int nd, int ord = 2, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline void calc_grad(arrvec_t v, arr_t a, ijk_t ijk, dijk_t dijk, typename std::enable_if<nd == 1>::type* = 0)
      {
        v[0](ijk) = grad<0, ord>(a, ijk[0], dijk[0]);
      }

      // 2D version
      template <int nd, int ord = 2, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline void calc_grad(arrvec_t v, arr_t a, ijk_t ijk, dijk_t dijk, typename std::enable_if<nd == 2>::type* = 0)
      {
        v[0](ijk) = grad<0, ord>(a, ijk[0], ijk[1], dijk[0]);
        v[1](ijk) = grad<1, ord>(a, ijk[1], ijk[0], dijk[1]);
      }

      // 3D version
      template <int nd, int ord = 2, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline void calc_grad(arrvec_t v, arr_t a, ijk_t ijk, dijk_t dijk, typename std::enable_if<nd == 3>::type* = 0)
      {
        v[0](ijk) = grad<0, ord>(a, ijk[0], ijk[1], ijk[2], dijk[0]);
        v[1](ijk) = grad<1, ord>(a, ijk[1], ijk[2], ijk[0], dijk[1]);
        v[2](ijk) = grad<2, ord>(a, ijk[2], ijk[0], ijk[1], dijk[2]);
      }
      
      // 2D version
      template <int nd, class arrvec_t, class arr_t, class ijk_t, class ijkm_t, class dijk_t>
      inline void calc_grad_cmpct(arrvec_t v, arr_t a, ijk_t ijk, ijkm_t ijkm, dijk_t dijk, typename std::enable_if<nd == 2>::type* = 0)
      {
        v[0](ijkm[0] + h, ijk[1]) = formulae::nabla::grad_cmpct<0>(a, ijkm[0], ijk[1], dijk[0]);
        v[1](ijk[0], ijkm[1] + h) = formulae::nabla::grad_cmpct<1>(a, ijkm[1], ijk[0], dijk[1]);
      }

      // 3D version
      template <int nd, class arrvec_t, class arr_t, class ijk_t, class ijkm_t, class dijk_t>
      inline void calc_grad_cmpct(arrvec_t v, arr_t a, ijk_t ijk, ijkm_t ijkm, dijk_t dijk, typename std::enable_if<nd == 3>::type* = 0)
      {
        v[0](ijkm[0] + h, ijk[1], ijk[2]) = formulae::nabla::grad_cmpct<0>(a, ijkm[0], ijk[1], ijk[2], dijk[0]);
        v[1](ijk[0], ijkm[1] + h, ijk[2]) = formulae::nabla::grad_cmpct<1>(a, ijkm[1], ijk[2], ijk[0], dijk[1]);
        v[2](ijk[0], ijk[1], ijkm[2] + h) = formulae::nabla::grad_cmpct<2>(a, ijkm[2], ijk[0], ijk[1], dijk[2]);
      }
      
      // divergence
      
      // 2D version
      template <int nd, int ord = 2, class arrvec_t, class ijk_t, class dijk_t>
      inline auto div(
	const arrvec_t &v, // vector field
	const ijk_t &ijk,
	const dijk_t dijk,
        typename std::enable_if<nd == 2>::type* = 0
      ) return_macro(,
          grad<0 BOOST_PP_COMMA() ord>(v[0], ijk[0], ijk[1], dijk[0])
        + grad<1 BOOST_PP_COMMA() ord>(v[1], ijk[1], ijk[0], dijk[1])
      )
      
      // 3D version
      template <int nd, int ord = 2, class arrvec_t, class ijk_t, class dijk_t>
      inline auto div(
	const arrvec_t &v, // vector field
	const ijk_t &ijk,
	const dijk_t dijk,
        typename std::enable_if<nd == 3>::type* = 0
      ) return_macro(,
          grad<0 BOOST_PP_COMMA() ord>(v[0], ijk[0], ijk[1], ijk[2], dijk[0])
        + grad<1 BOOST_PP_COMMA() ord>(v[1], ijk[1], ijk[2], ijk[0], dijk[1])
        + grad<2 BOOST_PP_COMMA() ord>(v[2], ijk[2], ijk[0], ijk[1], dijk[2])
      )
    } // namespace nabla_op
  } // namespace formulae
} // namespace libmpdataxx
