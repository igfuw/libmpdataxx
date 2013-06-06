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
      using namespace arakawa_c;

      // 1D
// TODO: psi -> psi/rho !!!
      template <opts_t opts, class arr_1d_t>
      inline auto beta_up(
        const arr_1d_t &psi,
        const arr_1d_t &psi_max, // from before the first iteration
        const arr_1d_t &C_corr, 
        const rng_t i  
      ) return_macro(,
        frac<opts>(
            max(max(max(psi_max(i), psi(i-1)), psi(i)), psi(i+1)) 
          - psi(i)
          ,// ----------------------------
// opts.pdf version
//            max(0, C_corr(i-h)) * psi(i-1) 
//          - min(0, C_corr(i+h)) * psi(i+1)

            max(0, C_corr(i-h)) * max(0, psi(i-1))
          - min(0, C_corr(i+h)) * max(0, psi(i+1))
          - max(0, C_corr(i+h)) * min(0, psi(i))
          + min(0, C_corr(i-h)) * min(0, psi(i))
        ) 
      ) 

// TODO: psi -> psi/rho !!!
      template <opts_t opts, class arr_1d_t>
      inline auto beta_dn(
        const arr_1d_t &psi, 
        const arr_1d_t &psi_min, // from before the first iteration
        const arr_1d_t &C_corr, 
        const rng_t i 
      ) return_macro(,
        frac<opts>(
            psi(i)
          - min(min(min(psi_min(i), psi(i-1)), psi(i)), psi(i+1)) 
          ,// --------------------------
// opts.pdf version
//            max(0, C_corr(i+h)) * psi(i) 
//          - min(0, C_corr(i-h)) * psi(i)
            max(0, C_corr(i+h)) * max(0, psi(i))
          - min(0, C_corr(i-h)) * max(0, psi(i))
          - max(0, C_corr(i-h)) * min(0, psi(i-1))
          + min(0, C_corr(i+h)) * min(0, psi(i+1))
        ) 
      ) 

    }; // namespace mpdata_fct
  }; // namespace formulae
}; // namespcae libmpdataxx
