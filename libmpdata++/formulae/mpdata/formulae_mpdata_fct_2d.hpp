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

      // TODO: add formula for c_mono<d>


// TODO: psi -> psi/rho !!!
      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_up(
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arr_2d_t &C_corr, 
        const rng_t i,  
        const rng_t j
      ) return_macro(,
        frac<opts>(
            max(max(max(max(max(psi_max(pi<d>(i, j)), 
                                   psi(pi<d>(i, j+1))),
              psi(pi<d>(i-1, j))), psi(pi<d>(i, j  ))), psi(pi<d>(i+1, j))),
                                   psi(pi<d>(i, j-1))
            )
          - psi(pi<d>(i, j))
          ,// ----------------------------
// opts.pdf version
//            max(0, C_corr(i-h)) * psi(i-1) 
//          - min(0, C_corr(i+h)) * psi(i+1)

            max(0, C_corr(pi<d>(i-h, j))) * max(0, psi(pi<d>(i-1, j))) 
          - min(0, C_corr(pi<d>(i+h, j))) * max(0, psi(pi<d>(i+1, j)))
          - max(0, C_corr(pi<d>(i+h, j))) * min(0, psi(pi<d>(i,   j)))
          + min(0, C_corr(pi<d>(i-h, j))) * min(0, psi(pi<d>(i,   j)))

          + max(0, C_corr(pi<d>(i, j-h))) * max(0, psi(pi<d>(i, j-1))) // TODO: double check indices
          - min(0, C_corr(pi<d>(i, j+h))) * max(0, psi(pi<d>(i, j+1)))
          - max(0, C_corr(pi<d>(i, j+h))) * min(0, psi(pi<d>(i, j  )))
          + min(0, C_corr(pi<d>(i, j-h))) * min(0, psi(pi<d>(i, j  )))
        ) 
      ) 

// TODO: psi -> psi/rho !!!
      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_dn(
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arr_2d_t &C_corr, 
        const rng_t i,
        const rng_t j 
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
