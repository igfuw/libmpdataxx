/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>

namespace libmpdataxx
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // infinite gauge version (for both signed and variable-sign fields)
      template <opts_t opts, class arr_1d_t>
      inline auto A(
        const arr_1d_t &psi, 
        const rng_t &i,
        typename std::enable_if<opt_set(opts, iga)>::type* = 0 // enabled if iga == true
      ) return_macro(,
	( 
	    psi(i+1)
	  - psi(i  )
	) / (// ------
	    1
	  + 1
	)
      )

      // single-sign signal version
      template <opts_t opts, class arr_1d_t>
      inline auto A(
        const arr_1d_t &psi, 
        const rng_t &i,
        typename std::enable_if<!opt_set(opts, iga) && opt_set(opts, pds)>::type* = 0 
      ) return_macro(,
	frac<opts>( 
	    psi(i+1)
	  - psi(i  )
	  ,// ------
	    psi(i+1)
	  + psi(i  )
	)
      )

      // variable-sign signal version (likely a good default)
      template <opts_t opts, class arr_1d_t>
      inline auto A(
        const arr_1d_t &psi, 
        const rng_t &i,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0 
      ) return_macro(,
	frac<opts>( 
	    abs(psi(i+1)) 
	  - abs(psi(i  ))
	  ,// -----------
	    abs(psi(i+1)) 
	  + abs(psi(i  ))
	) 
      )


      // 3rd order term (first term from eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998 (with G=1))
      template<opts_t opts, class arr_1d_t>
      inline auto HOT(
        const arr_1d_t &psi,
        const arr_1d_t &C,
        const rng_t &i,
        typename std::enable_if<!opt_set(opts, toa)>::type* = 0 
      ) -> decltype(0)
      { 
        return 0; 
      }

      template<class arr_1d_t>
      inline auto HOT_helper(
        const arr_1d_t &C,
        const rng_t &i
      ) return_macro(,
	(3 * C(i+h) * abs(C(i+h)) - 2 * pow(C(i+h), 3) - C(i+h)) / 3
      )

      template<opts_t opts, class arr_1d_t>
      inline auto HOT(
        const arr_1d_t &psi,
        const arr_1d_t &C,
        const rng_t &i,
        typename std::enable_if<opt_set(opts, toa) && opt_set(opts, pds) && !opt_set(opts, iga)>::type* = 0 
      ) return_macro(,
        HOT_helper(C, i) 
        * frac<opts>(
	  psi(i+2) - psi(i+1) - psi(i) + psi(i-1)
	  ,// 
	  psi(i+2) + psi(i+1) + psi(i) + psi(i-1)
	)
      )

      template<opts_t opts, class arr_1d_t>
      inline auto HOT(
        const arr_1d_t &psi,
        const arr_1d_t &C,
        const rng_t &i,
        typename std::enable_if<opt_set(opts, toa) && !opt_set(opts, pds) && !opt_set(opts, iga)>::type* = 0
      ) return_macro(,
        HOT_helper(C, i)
        * frac<opts>(
	  abs(psi(i+2)) - abs(psi(i+1)) - abs(psi(i)) + abs(psi(i-1))
	  ,// 
	  abs(psi(i+2)) + abs(psi(i+1)) + abs(psi(i)) + abs(psi(i-1))
	)
      )

      template<opts_t opts, class arr_1d_t>
      inline auto HOT(
        const arr_1d_t &psi,
        const arr_1d_t &C,
        const rng_t &i,
        typename std::enable_if<opt_set(opts, toa) && !opt_set(opts, pds) && opt_set(opts, iga)>::type* = 0 
      ) return_macro(,
        HOT_helper(C, i) 
        * (
	  psi(i+2) - psi(i+1) - psi(i) + psi(i-1)
	) / (
	  1 + 1 + 1 + 1
	)
      )

      //
      template<opts_t opts, class arr_1d_t>
      inline auto antidiff(
        const arr_1d_t &psi, 
        const rng_t &i, 
        const arr_1d_t &C
      ) return_macro(,
        // second-order terms
        abs(C(i+h)) 
        * (1 - abs(C(i+h))) 
        * A<opts>(psi, i) 
        // third-order terms
        + 
        HOT<opts>(psi, C, i) 
      )
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
