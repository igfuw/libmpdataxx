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

// TODO: psi -> psi/rho !!!
      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_up(
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
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
// TODO: opts.pds version
//            max(0, C_corr(i-h)) * psi(i-1) 
//          - min(0, C_corr(i+h)) * psi(i+1)

            pospart<opts>(C_corr[d-0](pi<d>(i-h, j))) * pospart<opts>(psi(pi<d>(i-1, j))) 
          - negpart<opts>(C_corr[d-0](pi<d>(i+h, j))) * pospart<opts>(psi(pi<d>(i+1, j)))
          - pospart<opts>(C_corr[d-0](pi<d>(i+h, j))) * negpart<opts>(psi(pi<d>(i,   j)))
          + negpart<opts>(C_corr[d-0](pi<d>(i-h, j))) * negpart<opts>(psi(pi<d>(i,   j)))

          + pospart<opts>(C_corr[d-1](pi<d>(i, j-h))) * pospart<opts>(psi(pi<d>(i, j-1))) // TODO: double check indices
          - negpart<opts>(C_corr[d-1](pi<d>(i, j+h))) * pospart<opts>(psi(pi<d>(i, j+1)))
          - pospart<opts>(C_corr[d-1](pi<d>(i, j+h))) * negpart<opts>(psi(pi<d>(i, j  )))
          + negpart<opts>(C_corr[d-1](pi<d>(i, j-h))) * negpart<opts>(psi(pi<d>(i, j  )))
        ) 
      ) 

// TODO: psi -> psi/rho !!!
      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_dn(
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,
        const rng_t j 
      ) return_macro(,
        frac<opts>(
            psi(pi<d>(i, j))
          - min(min(min(min(min(psi_min(pi<d>(i,j)), 
                                   psi(pi<d>(i, j+1))),
              psi(pi<d>(i-1, j))), psi(pi<d>(i, j  ))), psi(pi<d>(i+1, j))), // TODO: this is repeated with other dimension
                                   psi(pi<d>(i, j-1))
          )
          ,// --------------------------
// TODO: opts.pds version
//            max(0, C_corr(i+h)) * psi(i) 
//          - min(0, C_corr(i-h)) * psi(i)
            pospart<opts>(C_corr[d-0](pi<d>(i+h, j))) * pospart<opts>(psi(pi<d>(i,   j)))
          - negpart<opts>(C_corr[d-0](pi<d>(i-h, j))) * pospart<opts>(psi(pi<d>(i,   j)))
          - pospart<opts>(C_corr[d-0](pi<d>(i-h, j))) * negpart<opts>(psi(pi<d>(i-1, j)))
          + negpart<opts>(C_corr[d-0](pi<d>(i+h, j))) * negpart<opts>(psi(pi<d>(i+1, j)))

          + pospart<opts>(C_corr[d-1](pi<d>(i, j+h))) * pospart<opts>(psi(pi<d>(i,   j)))
          - negpart<opts>(C_corr[d-1](pi<d>(i, j-h))) * pospart<opts>(psi(pi<d>(i,   j)))
          - pospart<opts>(C_corr[d-1](pi<d>(i, j-h))) * negpart<opts>(psi(pi<d>(i, j-1)))
          + negpart<opts>(C_corr[d-1](pi<d>(i, j+h))) * negpart<opts>(psi(pi<d>(i, j+1)))
        ) 
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto C_mono(
        const arr_2d_t &psi,
        const arr_2d_t &psi_min, // from before the first iteration
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,
        const rng_t j
      ) return_macro(,
        C_corr[d]( pi<d>(i+h, j) ) * where(
          // if
          C_corr[d]( pi<d>(i+h, j) ) > 0,
          // then
          min(1, min(
            beta_dn<opts, d>(psi, psi_min, C_corr, i,     j),
            beta_up<opts, d>(psi, psi_max, C_corr, i + 1, j)
          )),
          // else
          min(1, min(
            beta_up<opts, d>(psi, psi_max, C_corr, i,     j),
            beta_dn<opts, d>(psi, psi_min, C_corr, i + 1, j)
          ))
        )
      ) 

    }; // namespace mpdata_fct
  }; // namespace formulae
}; // namespcae libmpdataxx
