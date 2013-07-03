/** @file
  * @copyright University of Warsaw
  * @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
  * @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
  * @brief Flux Corrected Transport formulae for MPDATA 
  *        (aka non-oscillatory, monotonic, sign-preserving option)
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
      /// \f$ \beta^{\uparrow}_{i} = \frac { \psi^{max}_{i}- \psi^{*}_{i} }
      /// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i-1/2}]^{+} \psi^{*}_{i-1} - 
      /// [u^{I}_{i+1/2}]^{-} \psi^{*}_{i+1} \right)  } \f$ \n
      /// eq.(19a) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
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
            pospart<opts>(C_corr(i-h)) * pospart<opts>(psi(i-1))
          - negpart<opts>(C_corr(i+h)) * pospart<opts>(psi(i+1))
          - pospart<opts>(C_corr(i+h)) * negpart<opts>(psi(i))
          + negpart<opts>(C_corr(i-h)) * negpart<opts>(psi(i))
        ) 
      ) 

// TODO: psi -> psi/rho !!!
      /// \f$ \beta^{\downarrow}_{i} = \frac { \psi^{*}_{i}- \psi^{min}_{i} }
      /// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i+1/2}]^{+} \psi^{*}_{i} - 
      /// [u^{I}_{i-1/2}]^{-} \psi^{*}_{i} \right)  } \f$ \n
      /// eq.(19b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
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
            pospart<opts>(C_corr(i+h)) * pospart<opts>(psi(i  ))
          - negpart<opts>(C_corr(i-h)) * pospart<opts>(psi(i  ))
          - pospart<opts>(C_corr(i-h)) * negpart<opts>(psi(i-1))
          + negpart<opts>(C_corr(i+h)) * negpart<opts>(psi(i+1))
        ) 
      ) 

      /// nonoscillatory antidiffusive velocity: \n
      /// \f$ U^{MON}_{i+1/2}=min(1,\beta ^{\downarrow}_i,\beta ^{\uparrow} _{i+1})[U_{i+1/2}]^{+} 
      /// + min(1,\beta^{\uparrow}_{i},\beta^{\downarrow}_{i+1/2})[u_{i+1/2}]^{-} \f$ \n
      /// where \f$ [\cdot]^{+}=max(\cdot,0) \f$ and \f$ [\cdot]^{-}=min(\cdot,0) \f$ \n
      /// eq.(18) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
      template <opts_t opts, class arr_1d_t>
      inline auto C_mono(
        const arr_1d_t &psi,
        const arr_1d_t &psi_min, // from before the first iteration
        const arr_1d_t &psi_max, // from before the first iteration
        const arr_1d_t &C_corr,
        const rng_t i
      ) return_macro(,
/* pds version TODO!!!
        C_corr( i+h ) * where(
          // if
          C_corr( i+h ) > 0,
          // then
          min(1, min(
            beta_dn<opts>(psi, psi_min, C_corr, i),
            beta_up<opts>(psi, psi_max, C_corr, i + 1)
          )), 
          // else
          min(1, min(
            beta_up<opts>(psi, psi_max, C_corr, i),
            beta_dn<opts>(psi, psi_min, C_corr, i + 1)
          ))  
        )  
*/
       C_corr(i+h) * where(
          // if
          C_corr(i+h) > 0,
          // then
          where(
            // if
            psi(i) > 0,
            // then
            min(1, min(
              beta_dn<opts>(psi, psi_min, C_corr, i    ), 
              beta_up<opts>(psi, psi_max, C_corr, i + 1)
            )), 
            // else
            min(1, min(
              beta_up<opts>(psi, psi_max, C_corr, i    ), 
              beta_dn<opts>(psi, psi_min, C_corr, i + 1)
            ))  
          ),  
          // else
          where(
            // if
            psi(i+1) > 0, // TODO: what if crossing zero?
            // then
            min(1, min(
              beta_up<opts>(psi, psi_max, C_corr, i    ), 
              beta_dn<opts>(psi, psi_min, C_corr, i + 1)
            )), 
            // else
            min(1, min(
              beta_dn<opts>(psi, psi_min, C_corr, i   ), 
              beta_up<opts>(psi, psi_max, C_corr, i + 1)
            ))  
          )   
        )   
      )
    }; // namespace mpdata_fct
  }; // namespace formulae
}; // namespcae libmpdataxx
