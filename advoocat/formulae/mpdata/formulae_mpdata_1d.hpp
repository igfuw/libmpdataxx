/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <advoocat/formulae/mpdata/formulae_mpdata_common.hpp>

namespace advoocat 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // 1D
      template <class arr_1d_t>
      inline auto A(
        const arr_1d_t &psi, 
        const rng_t &i 
      ) return_macro(,
        frac(
            abs(psi(i+1)) 
          - abs(psi(i  )),
          // ----------------------
            abs(psi(i+1)) 
          + abs(psi(i  ))
        ) 
      ) 

      template<class arr_1d_t>
      inline auto antidiff(
        const arr_1d_t &psi, 
        const rng_t &i, 
        const arr_1d_t &C
      ) return_macro(,
        abs(C(i+h)) 
        * (1 - abs(C(i+h))) 
        * A(psi, i) 
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae advoocat 
