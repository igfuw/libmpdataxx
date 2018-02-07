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

      template <class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const real_t dx
      )
      {
        return blitz::safeToReturn(
  	  (
	    x(i+1) - 
	    x(i-1)
	  ) / dx / 2.
        );
      }

      // 2D version
      template <int d, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const real_t dx
      )
      {
        return blitz::safeToReturn(
  	  (
  	    x(pi<d>(i+1, j)) - 
  	    x(pi<d>(i-1, j))
  	  ) / dx / 2.
        );
      }
      
      // 3D version
      template <int d, class arg_t, typename real_t>
      inline auto grad(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k,
	const real_t dx
      ) 
      {
        return blitz::safeToReturn(
  	  (
	    x(pi<d>(i+1, j, k)) - 
	    x(pi<d>(i-1, j, k))
	  ) / dx / 2.
        );
      }
      
      template <class arg_t, typename real_t>
      inline auto grad_cmpct(
	const arg_t &x,
	const rng_t &i,
	const real_t dx
      )
      {
        return blitz::safeToReturn(
	  (
	    x(i+1) - 
	    x(i)
	  ) / dx
        );
      }

      // 2D version
      template <int d, class arg_t, typename real_t>
      inline auto grad_cmpct(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const real_t dx
      ) 
      {
        return blitz::safeToReturn(
	  (
	    x(pi<d>(i+1, j)) - 
	    x(pi<d>(i  , j))
	  ) / dx
        );
      }
      
      // 3D version
      template <int d, class arg_t, typename real_t>
      inline auto grad_cmpct(
	const arg_t &x,
	const rng_t &i,
	const rng_t &j,
	const rng_t &k,
	const real_t dx
      )
      {
        return blitz::safeToReturn(
	  (
	    x(pi<d>(i+1, j, k)) - 
	    x(pi<d>(i, j, k))
	  ) / dx
        );
      }
      
      // helper function to calculate gradient components of a scalar field
      
      // 1D version
      template <int nd, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline void calc_grad(arrvec_t v, arr_t a, ijk_t ijk, dijk_t dijk, typename std::enable_if<nd == 1>::type* = 0)
      {
        v[0](ijk) = formulae::nabla::grad<0>(a, ijk[0], dijk[0]);
      }

      // 2D version
      template <int nd, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline void calc_grad(arrvec_t v, arr_t a, ijk_t ijk, dijk_t dijk, typename std::enable_if<nd == 2>::type* = 0)
      {
        v[0](ijk) = formulae::nabla::grad<0>(a, ijk[0], ijk[1], dijk[0]);
        v[1](ijk) = formulae::nabla::grad<1>(a, ijk[1], ijk[0], dijk[1]);
      }

      // 3D version
      template <int nd, class arrvec_t, class arr_t, class ijk_t, class dijk_t>
      inline void calc_grad(arrvec_t v, arr_t a, ijk_t ijk, dijk_t dijk, typename std::enable_if<nd == 3>::type* = 0)
      {
        v[0](ijk) = formulae::nabla::grad<0>(a, ijk[0], ijk[1], ijk[2], dijk[0]);
        v[1](ijk) = formulae::nabla::grad<1>(a, ijk[1], ijk[2], ijk[0], dijk[1]);
        v[2](ijk) = formulae::nabla::grad<2>(a, ijk[2], ijk[0], ijk[1], dijk[2]);
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
      template <int nd, class arrvec_t, class ijk_t, class dijk_t>
      inline auto div(
	const arrvec_t &v, // vector field
	const ijk_t &ijk,
	const dijk_t dijk,
        typename std::enable_if<nd == 2>::type* = 0
      ) 
      {
        return blitz::safeToReturn(
	  (v[0](ijk[0]+1, ijk[1]) - v[0](ijk[0]-1, ijk[1])) / dijk[0] / 2.
	  +
	  (v[1](ijk[0], ijk[1]+1) - v[1](ijk[0], ijk[1]-1)) / dijk[1] / 2.
        );
      }
      
      // 3D version
      template <int nd, class arrvec_t, class ijk_t, class dijk_t>
      inline auto div(
	const arrvec_t &v, // vector field
	const ijk_t &ijk,
	const dijk_t dijk,
        typename std::enable_if<nd == 3>::type* = 0
      ) 
      {
        return blitz::safeToReturn(
	  (v[0](ijk[0]+1, ijk[1], ijk[2]) - v[0](ijk[0]-1, ijk[1], ijk[2])) / dijk[0] / 2.
	  +
	  (v[1](ijk[0], ijk[1]+1, ijk[2]) - v[1](ijk[0], ijk[1]-1, ijk[2])) / dijk[1] / 2.
	  +
	  (v[2](ijk[0], ijk[1], ijk[2]+1) - v[2](ijk[0], ijk[1], ijk[2]-1)) / dijk[2] / 2.
        );
      }
    } // namespace nabla_op
  } // namespace formulae
} // namespace libmpdataxx
