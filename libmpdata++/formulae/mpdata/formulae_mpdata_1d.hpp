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
        where(opt_set(opts, sss),
          frac( // single-sign signal version
	      psi(i+1)
	    - psi(i  )
	    ,// ------
	      psi(i+1)
	    + psi(i  )
	  ),
	  frac( // variable-sign signal version (likely a good default)
	      abs(psi(i+1)) 
	    - abs(psi(i  ))
	    ,// -----------
	      abs(psi(i+1)) 
	    + abs(psi(i  ))
	  ) 
        )
      ) 

      template<opts_t opts, class arr_1d_t>
      inline auto antidiff(
        const arr_1d_t &psi, 
        const rng_t &i, 
        const arr_1d_t &C
      ) return_macro(,
        abs(C(i+h)) 
        * (1 - abs(C(i+h))) 
        * A<opts>(psi, i) 
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx
