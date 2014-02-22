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

      template<opts_t opts, class arr_1d_t>
      inline auto beta_up_nominator(
        const arr_1d_t &psi,
        const arr_1d_t &psi_max, 
        const arr_1d_t &G,
        const rng_t i
      ) return_macro(,
        (max(max(max(psi_max(i), psi(i-1)), psi(i)), psi(i+1)) - psi(i)) * formulae::G<opts>(G, i)
      )                                                                    //to make beta dimensionless 
                                                                           //when transporting mixing ratios with momentum
      template <opts_t opts, class arr_1d_t>
      inline auto beta_up( // for positive sign signal
        const arr_1d_t &psi,
        const arr_1d_t &psi_max, // from before the first iteration
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t i,
        typename std::enable_if<!opts::isset(opts, opts::abs) && !opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_up_nominator<opts>(psi, psi_max, G, i)
          ,// ----------------------------
            pospart<opts>(GC_corr(i-h)) * psi(i-1) 
          - negpart<opts>(GC_corr(i+h)) * psi(i+1)
        ) 
      ) 

      template <opts_t opts, class arr_1d_t>
      inline auto beta_up( //for variable sign signal
        const arr_1d_t &psi,
        const arr_1d_t &psi_max, // from before the first iteration
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_up_nominator<opts>(psi, psi_max, G, i)
          ,// ----------------------------
            pospart<opts>(GC_corr(i-h)) * pospart<opts>(psi(i-1)) // TODO: some parenthesis?
          - negpart<opts>(GC_corr(i+h)) * pospart<opts>(psi(i+1))
          - pospart<opts>(GC_corr(i+h)) * negpart<opts>(psi(i  ))
          + negpart<opts>(GC_corr(i-h)) * negpart<opts>(psi(i  ))
        ) 
      ) 

      template <opts_t opts, class arr_1d_t>
      inline auto beta_up( //infinite gauge option (* psi -> * 1)
        const arr_1d_t &psi,
        const arr_1d_t &psi_max, // from before the first iteration
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
        frac<opts>(
          beta_up_nominator<opts>(psi, psi_max, G, i)
          ,// ----------------------------
          pospart<opts>(GC_corr(i-h))   /* * 1 */
          - negpart<opts>(GC_corr(i+h)) /* * 1 */
          //+ blitz::epsilon(typename arr_1d_t::T_numtype(0))
        ) 
      ) 

      /// \f$ \beta^{\downarrow}_{i} = \frac { \psi^{*}_{i}- \psi^{min}_{i} }
      /// { \sum\limits_{I} \frac{\Delta t}{\Delta x^{I}} \left( [u^{I}_{i+1/2}]^{+} \psi^{*}_{i} - 
      /// [u^{I}_{i-1/2}]^{-} \psi^{*}_{i} \right)  } \f$ \n
      /// eq.(19b) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)

      template<opts_t opts, class arr_1d_t>
      inline auto beta_dn_nominator(
        const arr_1d_t &psi,
        const arr_1d_t &psi_min, 
        const arr_1d_t &G, 
        const rng_t i
      ) return_macro(,
        (psi(i) - min(min(min(psi_min(i), psi(i-1)), psi(i)), psi(i+1))) * formulae::G<opts>(G, i)
      )                                                                    //to make beta dimensionless 
                                                                           //when transporting mixing ratios with momentum
      template <opts_t opts, class arr_1d_t>
      inline auto beta_dn( //positive sign signal
        const arr_1d_t &psi, 
        const arr_1d_t &psi_min, // from before the first iteration
        const arr_1d_t &GC_corr, 
        const arr_1d_t &G, 
        const rng_t i, 
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_dn_nominator<opts>(psi, psi_min, G, i)
          ,// --------------------------
            pospart<opts>(GC_corr(i+h)) * psi(i) 
          - negpart<opts>(GC_corr(i-h)) * psi(i)
        ) 
      ) 

      template <opts_t opts, class arr_1d_t>
      inline auto beta_dn( // variable sign signal
        const arr_1d_t &psi, 
        const arr_1d_t &psi_min, // from before the first iteration
        const arr_1d_t &GC_corr, 
        const arr_1d_t &G, 
        const rng_t i, 
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_dn_nominator<opts>(psi, psi_min, G, i)
          ,// --------------------------
            pospart<opts>(GC_corr(i+h)) * pospart<opts>(psi(i  ))
          - negpart<opts>(GC_corr(i-h)) * pospart<opts>(psi(i  ))
          - pospart<opts>(GC_corr(i-h)) * negpart<opts>(psi(i-1))
          + negpart<opts>(GC_corr(i+h)) * negpart<opts>(psi(i+1))
        ) 
      ) 

      template <opts_t opts, class arr_1d_t>
      inline auto beta_dn( // infinite gauge option
        const arr_1d_t &psi, 
        const arr_1d_t &psi_min, // from before the first iteration
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0 // enabled if iga == true
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs are mutually exclusive");
        ,
        frac<opts>(
          beta_dn_nominator<opts>(psi, psi_min, G, i)
          ,// --------------------------
            pospart<opts>(GC_corr(i+h)) /* * 1 */
          - negpart<opts>(GC_corr(i-h)) /* * 1 */
          //+ blitz::epsilon(typename arr_1d_t::T_numtype(0))
        ) 
      ) 

      /// nonoscillatory antidiffusive velocity: \n
      /// \f$ U^{MON}_{i+1/2}=min(1,\beta ^{\downarrow}_i,\beta ^{\uparrow} _{i+1})[U_{i+1/2}]^{+} 
      /// + min(1,\beta^{\uparrow}_{i},\beta^{\downarrow}_{i+1/2})[u_{i+1/2}]^{-} \f$ \n
      /// where \f$ [\cdot]^{+}=max(\cdot,0) \f$ and \f$ [\cdot]^{-}=min(\cdot,0) \f$ \n
      /// eq.(18) in Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
      template <opts_t opts, class arr_1d_t>
      inline auto GC_mono(  // for variable sign signal and no inf. gauge
        const arr_1d_t &psi,
        const arr_1d_t &psi_min, // from before the first iteration
        const arr_1d_t &psi_max, // from before the first iteration
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0 // enabled if iga == false
      ) return_macro(,
       GC_corr(i+h) * where(
          // if
          GC_corr(i+h) > 0,
          // then
          where(
            // if
            psi(i) > 0,
            // then
            min(1, min(
              beta_dn<opts>(psi, psi_min, GC_corr, G, i    ), 
              beta_up<opts>(psi, psi_max, GC_corr, G, i + 1)
            )), 
            // else
            min(1, min(
              beta_up<opts>(psi, psi_max, GC_corr, G, i    ), 
              beta_dn<opts>(psi, psi_min, GC_corr, G, i + 1)
            ))  
          ),  
          // else
          where(
            // if
            psi(i+1) > 0, // TODO: what if crossing zero?
            // then
            min(1, min(
              beta_up<opts>(psi, psi_max, GC_corr, G, i    ), 
              beta_dn<opts>(psi, psi_min, GC_corr, G, i + 1)
            )), 
            // else
            min(1, min(
              beta_dn<opts>(psi, psi_min, GC_corr, G, i   ), 
              beta_up<opts>(psi, psi_max, GC_corr, G, i + 1)
            ))  
          )   
        )   
      )

      template <opts_t opts, class arr_1d_t>
      inline auto GC_mono( // for inf. gauge option or positive sign signal
        const arr_1d_t &psi,
        const arr_1d_t &psi_min, // from before the first iteration
        const arr_1d_t &psi_max, // from before the first iteration
        const arr_1d_t &GC_corr,
        const arr_1d_t &G,
        const rng_t i,
        typename std::enable_if<opts::isset(opts, opts::iga) || !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        GC_corr( i+h ) * where(
          // if
          GC_corr( i+h ) > 0,
          // then
          min(1, min(
            beta_dn<opts>(psi, psi_min, GC_corr, G, i),
            beta_up<opts>(psi, psi_max, GC_corr, G, i + 1)
          )), 
          // else
          min(1, min(
            beta_up<opts>(psi, psi_max, GC_corr, G, i),
            beta_dn<opts>(psi, psi_min, GC_corr, G, i + 1)
          ))  
        )
      )  

    }; // namespace mpdata_fct
  }; // namespace formulae
}; // namespcae libmpdataxx
