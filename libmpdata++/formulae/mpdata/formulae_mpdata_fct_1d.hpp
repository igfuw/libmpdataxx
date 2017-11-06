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

      /// \f$ \beta^{\uparrow}_{i} = \frac { \psi^{max}_{i}- \psi^{*}_{i} }
      /// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i-1/2}]^{+} \psi^{*}_{i-1} - 
      /// [u^{I}_{i+1/2}]^{-} \psi^{*}_{i+1} \right)  } \f$ \n
      /// eq.(19a) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

      template<opts_t opts, class arr_1d_t, class ix_t>
      forceinline_macro auto beta_up_nominator(
        const arr_1d_t &psi,
        const arr_1d_t &psi_max, 
        const arr_1d_t &G,
        const ix_t &i
      )
      {
        return return_helper<ix_t>(
          (max<ix_t>(psi_max(i), psi(i-1), psi(i), psi(i+1)) - psi(i)) * formulae::G<opts>(G, i)
        );
      }                                                                    //to make beta dimensionless 
                                                                           //when transporting mixing ratios with momentum
      template <opts_t opts, class arr_1d_t, class flx_t>
      forceinline_macro void beta_up( // for positive sign signal
        arr_1d_t &b,
        const arr_1d_t &psi,
        const arr_1d_t &psi_max, // from before the first iteration
        const flx_t &flx,
        const arr_1d_t &G,
        const rng_t &ir
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          b(i) = 
          fct_frac<ix_t>(
            beta_up_nominator<opts>(psi, psi_max, G, i)
            , // ----------------------------
              pospart<opts, ix_t>(flx[0](i-h))
            - negpart<opts, ix_t>(flx[0](i+h))
          );
        }
      } 

      /// \f$ \beta^{\downarrow}_{i} = \frac { \psi^{*}_{i}- \psi^{min}_{i} }
      /// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i+1/2}]^{+} \psi^{*}_{i} - 
      /// [u^{I}_{i-1/2}]^{-} \psi^{*}_{i} \right)  } \f$ \n
      /// eq.(19b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

      template<opts_t opts, class arr_1d_t, class ix_t>
      forceinline_macro auto beta_dn_nominator(
        const arr_1d_t &psi,
        const arr_1d_t &psi_min, 
        const arr_1d_t &G, 
        const ix_t &i
      )
      {
        return return_helper<ix_t>(
          (psi(i) - min<ix_t>(psi_min(i), psi(i-1), psi(i), psi(i+1))) * formulae::G<opts>(G, i)
        );
      }                                                                      //to make beta dimensionless 
                                                                             //when transporting mixing ratios with momentum
      template <opts_t opts, class arr_1d_t, class flx_t>
      forceinline_macro void beta_dn( //positive sign signal
        arr_1d_t &b, 
        const arr_1d_t &psi, 
        const arr_1d_t &psi_min, // from before the first iteration
        const flx_t &flx, 
        const arr_1d_t &G, 
        const rng_t &ir
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          b(i) = 
          fct_frac<ix_t>(
            beta_dn_nominator<opts>(psi, psi_min, G, i)
            , // --------------------------
              pospart<opts, ix_t>(flx[0](i+h))
            - negpart<opts, ix_t>(flx[0](i-h))
          );
        }
      } 

      /// nonoscillatory antidiffusive velocity: \n
      /// \f$ U^{MON}_{i+1/2}=min(1,\beta ^{\downarrow}_i,\beta ^{\uparrow} _{i+1})[U_{i+1/2}]^{+} 
      /// + min(1,\beta^{\uparrow}_{i},\beta^{\downarrow}_{i+1/2})[u_{i+1/2}]^{-} \f$ \n
      /// where<ix_t> \f$ [\cdot]^{+}=max(\cdot,0) \f$ and \f$ [\cdot]^{-}=min(\cdot,0) \f$ \n
      /// eq.(18) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
      template <opts_t opts, class arr_1d_t>
      forceinline_macro void GC_mono(  // for variable sign signal and no inf. gauge
        arr_1d_t &GC_m,
        const arr_1d_t &psi,
        const arr_1d_t &beta_up,
        const arr_1d_t &beta_dn,
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t &ir,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0 // enabled if iga == false
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          GC_m(i+h) =
          GC_corr(i+h) * where<ix_t>(
             // if
             GC_corr(i+h) > 0,
             // then
             where<ix_t>(
               // if
               psi(i) > 0,
               // then
               min<ix_t>(1,
                 beta_dn(i    ), 
                 beta_up(i + 1)
               ), 
               // else
               min<ix_t>(1,
                 beta_up(i    ), 
                 beta_dn(i + 1)
               )  
             ),  
             // else
             where<ix_t>(
               // if
               psi(i+1) > 0,
               // then
               min<ix_t>(1,
                 beta_up(i    ), 
                 beta_dn(i + 1)
               ), 
               // else
               min<ix_t>(1,
                 beta_dn(i   ), 
                 beta_up(i + 1)
               )  
             )   
          );
        } 
      }

      template <opts_t opts, class arr_1d_t>
      forceinline_macro void GC_mono( // for inf. gauge option or positive sign signal
        arr_1d_t &GC_m,
        const arr_1d_t &psi,
        const arr_1d_t &beta_up,
        const arr_1d_t &beta_dn,
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t &ir,
        typename std::enable_if<opts::isset(opts, opts::iga) || !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          GC_m(i+h) = 
          GC_corr(i+h) * where<ix_t>(
            // if
            GC_corr(i+h) > 0,
            // then
            min<ix_t>(1,
              beta_dn(i),
              beta_up(i + 1)
            ), 
            // else
            min<ix_t>(1,
              beta_up(i),
              beta_dn(i + 1)
            )  
          );
        }
      }  
    } // namespace mpdata_fct
  } // namespace formulae
} // namespcae libmpdataxx
