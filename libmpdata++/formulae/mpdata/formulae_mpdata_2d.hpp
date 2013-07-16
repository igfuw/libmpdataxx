/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

// TODO: iga

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // single-sign scalar version
      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<opt_set(opts, pds)>::type* = 0 // enabled if pds == true
      ) return_macro(,
	frac<opts>(
	    psi(pi<d>(i+1, j)) 
	  - psi(pi<d>(i,   j))
	  ,// ----------------
	    psi(pi<d>(i+1, j))
	  + psi(pi<d>(i,   j))
	)
      )

      // variable-sign scalar version
      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, pds)>::type* = 0 // enabled if pds == false
      ) return_macro(,
	frac<opts>(
	    abs(psi(pi<d>(i+1, j))) 
	  - abs(psi(pi<d>(i,   j)))
	  ,// ---------------------
	    abs(psi(pi<d>(i+1, j))) 
	  + abs(psi(pi<d>(i,   j)))
	) 
      ) 

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B(
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<opt_set(opts, pds)>::type* = 0 // enabled if pds == true
      ) return_macro(,
	frac<opts>( 
	    psi(pi<d>(i+1, j+1))
	  + psi(pi<d>(i,   j+1)) 
	  - psi(pi<d>(i+1, j-1)) 
	  - psi(pi<d>(i,   j-1))
	  ,// ------------------
	    psi(pi<d>(i+1, j+1)) 
	  + psi(pi<d>(i,   j+1)) 
	  + psi(pi<d>(i+1, j-1)) 
	  + psi(pi<d>(i,   j-1))
	) / 2
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B(
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opt_set(opts, pds)>::type* = 0 // enabled if pds == false
      ) return_macro(,
	frac<opts>( 
	    abs(psi(pi<d>(i+1, j+1))) 
	  + abs(psi(pi<d>(i,   j+1))) 
	  - abs(psi(pi<d>(i+1, j-1))) 
	  - abs(psi(pi<d>(i,   j-1)))
	  ,// -----------------------
	    abs(psi(pi<d>(i+1, j+1))) 
	  + abs(psi(pi<d>(i,   j+1))) 
	  + abs(psi(pi<d>(i+1, j-1))) 
	  + abs(psi(pi<d>(i,   j-1)))
	) / 2
      )

      template<int d, class arr_2d_t>
      inline auto C_bar(
        const arr_2d_t &C, 
        const rng_t &i, 
        const rng_t &j
      ) return_macro(,
        (
          C(pi<d>(i+1, j+h)) + 
          C(pi<d>(i,   j+h)) +
          C(pi<d>(i+1, j-h)) + 
          C(pi<d>(i,   j-h)) 
        ) / 4
      )

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto antidiff(
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        const arrvec_t<arr_2d_t> &C
      ) return_macro(,
        abs(C[dim](pi<dim>(i+h, j))) 
        * (1 - abs(C[dim](pi<dim>(i+h, j)))) 
        * A<opts BOOST_PP_COMMA() dim>(psi, i, j) 
        - C[dim](pi<dim>(i+h, j)) 
        * C_bar<dim>(C[dim-1], i, j)
        * B<opts BOOST_PP_COMMA() dim>(psi, i, j)
        // TODO: toa
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx 
