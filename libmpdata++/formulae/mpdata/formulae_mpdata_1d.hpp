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
      // single-sign signal version
      template <opts_t opts, class arr_1d_t>
      inline auto A(
        const arr_1d_t &psi, 
        const rng_t &i,
        typename std::enable_if<opt_set(opts, sss)>::type* = 0 // enabled if sss == true
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
        typename std::enable_if<!opt_set(opts, sss)>::type* = 0 // enabled if sss = false
      ) return_macro(,
	frac<opts>( 
	    abs(psi(i+1)) 
	  - abs(psi(i  ))
	  ,// -----------
	    abs(psi(i+1)) 
	  + abs(psi(i  ))
	) 
      )

// TODO: move toa formulae to a separate file
// TODO: rename A, B nd xx_helper to 1/psi dpsi/dt etc


      // 3rd order term (first term from eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998 (with G=1))
      template<opts_t opts, class arr_1d_t>
      inline auto f_3rd_xx(
        const arr_1d_t &psi,
        const arr_1d_t &C,
        const rng_t &i,
        typename std::enable_if<opt_set(opts, sss)>::type* = 0 // enabled if sss == true
      ) return_macro(,
	(3 * C(i+h) * abs(C(i+h)) - 2 * pow(C(i+h), 3) - C(i+h)) 
	/ (typename arr_1d_t::T_numtype(3)) // TODO: check if this cast is needed
        * // TODO: code duplication above :*
        frac<opts>(
	  psi(i+2) - psi(i+1) - psi(i) + psi(i-1)
	  ,// 
	  psi(i+2) + psi(i+1) + psi(i) + psi(i-1)
	)
    )

      template<opts_t opts, class arr_1d_t>
      inline auto f_3rd_xx(
        const arr_1d_t &psi,
        const arr_1d_t &C,
        const rng_t &i,
        typename std::enable_if<!opt_set(opts, sss)>::type* = 0 // enabled if sss == false
      ) return_macro(,
	(3 * C(i+h) * abs(C(i+h)) - 2 * pow(C(i+h), 3) - C(i+h)) 
	/ (typename arr_1d_t::T_numtype(3)) // TODO: check if this cast is needed
        * // TODO: code duplication above :*
        frac<opts>(
	  abs(psi(i+2)) - abs(psi(i+1)) - abs(psi(i)) + abs(psi(i-1))
	  ,// 
	  abs(psi(i+2)) + abs(psi(i+1)) + abs(psi(i)) + abs(psi(i-1))
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
        where(opt_set(opts, toa), f_3rd_xx<opts>(psi, C, i), 0) // TODO: rename to HOT? // TODO: where -> enable_if
      )
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
