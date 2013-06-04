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
      // 1D
      template <opts_t opts, class arr_1d_t>
      inline auto A(
        const arr_1d_t &psi, 
        const rng_t &i 
      ) return_macro(,
        where(
          // if
          opt_set(opts, sss),
          // then
          frac( // single-sign signal version
	      psi(i+1)
	    - psi(i  )
	    ,// ------
	      psi(i+1)
	    + psi(i  )
	  ),
          // else
	  frac( // variable-sign signal version (likely a good default)
	      abs(psi(i+1)) 
	    - abs(psi(i  ))
	    ,// -----------
	      abs(psi(i+1)) 
	    + abs(psi(i  ))
	  ) 
        )
      ) 

      // 3rd order term (first term from eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998 (with G=1))
      template<opts_t opts, class arr_1d_t>
      inline auto f_3rd_xx(
        const arr_1d_t &psi,
        const arr_1d_t &C,
        const rng_t &i
      ) return_macro(,
	(3 * C(i+h) * abs(C(i+h)) - 2 * pow(C(i+h), 3) - C(i+h)) 
        / (typename arr_1d_t::T_numtype(3)) // TODO: chec if this cast is needed
	* where(
          // if
          opt_set(opts, sss),
          // then
          frac(
	    psi(i+2) - psi(i+1) - psi(i) + psi(i-1)
	    ,// 
	    psi(i+2) + psi(i+1) + psi(i) + psi(i-1)
	  ),
          // else
          frac(
	    abs(psi(i+2)) - abs(psi(i+1)) - abs(psi(i)) + abs(psi(i-1))
	    ,// 
	    abs(psi(i+2)) + abs(psi(i+1)) + abs(psi(i)) + abs(psi(i-1))
	  )
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
        + where(
          // if
          opt_set(opts, toa),
          // then 
          f_3rd_xx<opts>(psi, C, i),
          // else
          0
        )
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
