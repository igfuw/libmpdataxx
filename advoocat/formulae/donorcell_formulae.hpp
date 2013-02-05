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
    namespace donorcell
    {
      using namespace arakawa_c;
      using idxperm::pi;

      template<class T1, class T2, class T3> 
      inline auto F(
	const T1 &psi_l, const T2 &psi_r, const T3 &C
      ) return_macro(,
	(
	  (C + abs(C)) * psi_l + 
	  (C - abs(C)) * psi_r
	) / 2
      )

      template <class arr_1d_t>
      inline auto donorcell( 
	const arr_1d_t &psi, 
	const arr_1d_t &C, 
	const rng_t &i
      ) return_macro(,
	F(
	  psi(i  ), 
	  psi(i+1), 
	    C(i+h)
	) -
	F(
	  psi(i-1), 
	  psi(i  ), 
	    C(i-h)
	)
      )

      template<int d, class arr_2d_t>  
      inline auto donorcell( 
	const arr_2d_t &psi, 
	const arr_2d_t &C, 
	const rng_t &i, 
	const rng_t &j
      ) return_macro(,
	F(
	  psi(pi<d>(i,   j)), 
	  psi(pi<d>(i+1, j)), 
	    C(pi<d>(i+h, j))
	) -
	F(
	  psi(pi<d>(i-1, j)), 
	  psi(pi<d>(i,   j)), 
	    C(pi<d>(i-h, j))
	)
      )

      template<int d, class arr_3d_t>  
      inline auto donorcell( 
	const arr_3d_t &psi, 
	const arr_3d_t &C, 
	const rng_t &i, 
	const rng_t &j,
	const rng_t &k
      ) return_macro(,
	F(
	  psi(pi<d>(i,   j, k)), 
	  psi(pi<d>(i+1, j, k)), 
	    C(pi<d>(i+h, j, k))
	) -
	F(
	  psi(pi<d>(i-1, j, k)), 
	  psi(pi<d>(i,   j, k)), 
	    C(pi<d>(i-h, j, k))
	)
      )

      template <class arr_1d_t>
      void op_1d(
	const arrvec_t<arr_1d_t> &psi, 
	const int n,
	const arr_1d_t &C, 
	const rng_t &i
      ) { 
	psi[n+1](i) = psi[n](i)
	  - donorcell(psi[n], C, i);
      }

      template <class arr_2d_t>
      void op_2d(
	const arrvec_t<arr_2d_t> &psi, const int n,
	const arrvec_t<arr_2d_t> &C, 
	const rng_t &i, const rng_t &j
      ) { 
	psi[n+1](i,j) = psi[n](i,j)
	  - donorcell<0>(psi[n], C[0], i, j)
	  - donorcell<1>(psi[n], C[1], j, i); 
      }

      template <class arr_3d_t>
      void op_3d(
	const arrvec_t<arr_3d_t> &psi, const int n,
	const arrvec_t<arr_3d_t> &C, 
	const rng_t &i, const rng_t &j, const rng_t &k
      ) { 
	psi[n+1](i,j) = psi[n](i,j)
	  - donorcell<0>(psi[n], C[0], i, j, k)
	  - donorcell<1>(psi[n], C[1], j, k, i)
	  - donorcell<2>(psi[n], C[2], k, i, j); 
      }
    }; // namespace donorcell 
  }; // namespace formulae
}; // namespace advoocat
