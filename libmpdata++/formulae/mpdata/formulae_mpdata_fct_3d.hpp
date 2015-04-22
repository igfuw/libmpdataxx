/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      using namespace arakawa_c;

      //see Smolarkiewicz & Grabowski 1990 (J.Comp.Phys.,86,355-375)
      template <opts_t opts, class arr_3d_t>
      inline auto beta_up_nominator(
        const arr_3d_t &psi,
        const arr_3d_t &psi_max,
        const arr_3d_t &G,
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k
      ) return_macro(, 
        (
          max(max(max(max(max(max(max(
                      psi_max(i, j, k),
                      psi(i, j, k)),
                      psi(i+1, j, k)),
                      psi(i-1, j, k)),
                      psi(i, j+1, k)),
                      psi(i, j-1, k)),
                      psi(i, j, k+1)),
                      psi(i, j, k-1)
          )
        - psi(i, j, k)
        ) * formulae::G<opts BOOST_PP_COMMA() 0>(G, i, j, k) //to make beta up dimensionless when transporting mixing ratios with momentum
      ) 

      template <opts_t opts, class arr_3d_t, class flx_t>
      inline auto beta_up(
        const arr_3d_t &psi,
        const arr_3d_t &psi_max, // from before the first iteration
        const flx_t &flx,
        const arr_3d_t &G,
        const rng_t &i, 
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
        fct_frac(
          beta_up_nominator<opts>(psi, psi_max, G, i, j, k)
          , //-------------------------------------------------------------------------------------
          ( pospart<opts>(flx[0](i-h, j, k))
          - negpart<opts>(flx[0](i+h, j, k)) )  // additional parenthesis so that we first sum
          +                                     // fluxes in separate dimensions
          ( pospart<opts>(flx[1](i, j-h, k))    // could be important for accuracy if one of them
          - negpart<opts>(flx[1](i, j+h, k)) )  // is of different magnitude than the other
          +                                     // fluxes in separate dimensions
          ( pospart<opts>(flx[2](i, j, k-h))
          - negpart<opts>(flx[2](i, j, k+h)) )
        )
      )

      template <opts_t opts, class arr_3d_t>
      inline auto beta_dn_nominator(
        const arr_3d_t &psi,
        const arr_3d_t &psi_min,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
        (
          psi(i, j, k)
          - min(min(min(min(min(min(min(
                    psi_min(i, j, k),
                    psi(i, j+1, k)),
                    psi(i-1, j, k)),
                    psi(i, j, k  )),
                    psi(i+1, j, k)),
                    psi(i, j, k+1)),
                    psi(i, j, k-1)),
                    psi(i, j-1, k)
          )
        ) * formulae::G<opts BOOST_PP_COMMA() 0>(G, i, j, k) //to make beta up dimensionless when transporting mixing ratios with momentum
      )
      
      template <opts_t opts, class arr_3d_t, class flx_t>
      inline auto beta_dn(
        const arr_3d_t &psi,
        const arr_3d_t &psi_min, // from before the first iteration
        const flx_t &flx,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
        fct_frac(
          beta_dn_nominator<opts>(psi, psi_min, G, i, j, k)
	  , //-----------------------------------------------------------------------------------
          ( pospart<opts>(flx[0](i+h, j, k))
          - negpart<opts>(flx[0](i-h, j, k)) )  //see note in beta up
          +
          ( pospart<opts>(flx[1](i, j+h, k))
          - negpart<opts>(flx[1](i, j-h, k)) )
          +
          ( pospart<opts>(flx[2](i, j, k+h))
          - negpart<opts>(flx[2](i, j, k-h)) )
        )
      )

      template <opts_t opts, int d, class arr_3d_t>
      inline auto GC_mono( //for variable-sign signal and no infinite gauge option
        const arr_3d_t &psi,
        const arr_3d_t &beta_up,
        const arr_3d_t &beta_dn,
        const arrvec_t<arr_3d_t> &GC_corr,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        GC_corr[d]( pi<d>(i+h, j, k) ) * where(
          // if
          GC_corr[d]( pi<d>(i+h, j, k) ) > 0,
          // then
          where(
            // if
            psi(pi<d>(i, j, k)) > 0,
            // then
	    min(1, min(
              beta_dn(pi<d>(i,     j, k)),
              beta_up(pi<d>(i + 1, j, k))
	    )),
            // else
	    min(1, min(
	      beta_up(pi<d>(i,     j, k)),
	      beta_dn(pi<d>(i + 1, j, k))
	    ))
          ),
          // else
          where(
            // if
            psi(pi<d>(i+1, j, k)) > 0,
            // then
	    min(1, min(
	      beta_up(pi<d>(i,     j, k)),
	      beta_dn(pi<d>(i + 1, j, k))
	    )),
            // else
	    min(1, min(
	      beta_dn(pi<d>(i,     j, k)),
	      beta_up(pi<d>(i + 1, j, k))
	    ))
          )
        )
      )

      template <opts_t opts, int d, class arr_3d_t>
      inline auto GC_mono( //for infinite gauge option or positive-sign signal
        const arr_3d_t &psi,
        const arr_3d_t &beta_up,
        const arr_3d_t &beta_dn,
        const arrvec_t<arr_3d_t> &GC_corr,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga) || !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        GC_corr[d]( pi<d>(i+h, j, k) ) * where(
          // if
          GC_corr[d]( pi<d>(i+h, j, k) ) > 0,
          // then
	  min(1, min(
	    beta_dn(pi<d>(i,     j, k)),
	    beta_up(pi<d>(i + 1, j, k))
	  )),
          // else
	  min(1, min(
	    beta_up(pi<d>(i,     j, k)),
	    beta_dn(pi<d>(i + 1, j, k))
	  ))
        )
      )
    }; // namespace mpdata_fct
  }; // namespace formulae
}; // namespcae libmpdataxx
