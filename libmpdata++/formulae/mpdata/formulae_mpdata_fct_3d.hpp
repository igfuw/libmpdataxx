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
      template <opts_t opts, class arr_3d_t, class ix_t>
      forceinline_macro auto beta_up_nominator(
        const arr_3d_t &psi,
        const arr_3d_t &psi_max,
        const arr_3d_t &G,
        const ix_t &i, 
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            max<ix_t>(psi_max(i, j, k),
                      psi(i, j, k),
                      psi(i+1, j, k),
                      psi(i-1, j, k),
                      psi(i, j+1, k),
                      psi(i, j-1, k),
                      psi(i, j, k+1),
                      psi(i, j, k-1)
            )
          - psi(i, j, k)
          ) * formulae::G<opts, 0>(G, i, j, k) //to make beta up dimensionless when transporting mixing ratios with momentum
        );
      } 

      template <opts_t opts, class arr_3d_t, class flx_t>
      forceinline_macro void beta_up(
        arr_3d_t &b,
        const arr_3d_t &psi,
        const arr_3d_t &psi_max, // from before the first iteration
        const flx_t &flx,
        const arr_3d_t &G,
        const rng_t &ir, 
        const rng_t &jr,
        const rng_t &kr
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            for (int k = kr.first(); k <= kr.last(); ++k)
            {
              b(i, j, k) = 
              fct_frac<ix_t>(
                beta_up_nominator<opts>(psi, psi_max, G, i, j, k)
                , //-------------------------------------------------------------------------------------
                ( pospart<opts, ix_t>(flx[0](i-h, j, k))
                - negpart<opts, ix_t>(flx[0](i+h, j, k)) )  // additional parenthesis so that we first sum
                +                                           // fluxes in separate dimensions
                ( pospart<opts, ix_t>(flx[1](i, j-h, k))    // could be important for accuracy if one of them
                - negpart<opts, ix_t>(flx[1](i, j+h, k)) )  // is of different magnitude than the other
                +                                           // fluxes in separate dimensions
                ( pospart<opts, ix_t>(flx[2](i, j, k-h))
                - negpart<opts, ix_t>(flx[2](i, j, k+h)) )
              );
            }
          }
        }
      }

      template <opts_t opts, class arr_3d_t, class ix_t>
      forceinline_macro auto beta_dn_nominator(
        const arr_3d_t &psi,
        const arr_3d_t &psi_min,
        const arr_3d_t &G,
        const ix_t &i,
        const ix_t &j,
        const ix_t &k
      )
      {
        return return_helper<ix_t>(
          (
            psi(i, j, k)
            - min<ix_t>(psi_min(i, j, k),
                        psi(i, j+1, k),
                        psi(i-1, j, k),
                        psi(i, j, k  ),
                        psi(i+1, j, k),
                        psi(i, j, k+1),
                        psi(i, j, k-1),
                        psi(i, j-1, k)
            )
          ) * formulae::G<opts, 0>(G, i, j, k) //to make beta up dimensionless when transporting mixing ratios with momentum
        );
      }
      
      template <opts_t opts, class arr_3d_t, class flx_t>
      forceinline_macro void beta_dn(
        arr_3d_t &b,
        const arr_3d_t &psi,
        const arr_3d_t &psi_min, // from before the first iteration
        const flx_t &flx,
        const arr_3d_t &G,
        const rng_t &ir,
        const rng_t &jr,
        const rng_t &kr
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            for (int k = kr.first(); k <= kr.last(); ++k)
            {
              b(i, j, k) = 
              fct_frac<ix_t>(
                beta_dn_nominator<opts>(psi, psi_min, G, i, j, k)
                , //-----------------------------------------------------------------------------------
                ( pospart<opts, ix_t>(flx[0](i+h, j, k))
                - negpart<opts, ix_t>(flx[0](i-h, j, k)) )  //see note in beta up
                +
                ( pospart<opts, ix_t>(flx[1](i, j+h, k))
                - negpart<opts, ix_t>(flx[1](i, j-h, k)) )
                +
                ( pospart<opts, ix_t>(flx[2](i, j, k+h))
                - negpart<opts, ix_t>(flx[2](i, j, k-h)) )
              );
            }
          }
        }
      }

      template <opts_t opts, int d, class arr_3d_t>
      forceinline_macro void GC_mono( //for variable-sign signal and no infinite gauge option
        arrvec_t<arr_3d_t> &GC_m,
        const arr_3d_t &psi,
        const arr_3d_t &beta_up,
        const arr_3d_t &beta_dn,
        const arrvec_t<arr_3d_t> &GC_corr,
        const arr_3d_t &G,
        const rng_t &ir,
        const rng_t &jr,
        const rng_t &kr,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            for (int k = kr.first(); k <= kr.last(); ++k)
            {
              GC_m[d]( pi<d>(i+h, j, k) ) =
              GC_corr[d]( pi<d>(i+h, j, k) ) * where<ix_t>(
                // if
                GC_corr[d]( pi<d>(i+h, j, k) ) > 0,
                // then
                where<ix_t>(
                  // if
                  psi(pi<d>(i, j, k)) > 0,
                  // then
                  min<ix_t>(1,
                    beta_dn(pi<d>(i,     j, k)),
                    beta_up(pi<d>(i + 1, j, k))
                  ),
                  // else
                  min<ix_t>(1,
                    beta_up(pi<d>(i,     j, k)),
                    beta_dn(pi<d>(i + 1, j, k))
                  )
                ),
                // else
                where<ix_t>(
                  // if
                  psi(pi<d>(i+1, j, k)) > 0,
                  // then
                  min<ix_t>(1,
                    beta_up(pi<d>(i,     j, k)),
                    beta_dn(pi<d>(i + 1, j, k))
                  ),
                  // else
                  min<ix_t>(1,
                    beta_dn(pi<d>(i,     j, k)),
                    beta_up(pi<d>(i + 1, j, k))
                  )
                )
              );
            }
          }
        }
      }

      template <opts_t opts, int d, class arr_3d_t>
      forceinline_macro void GC_mono( //for infinite gauge option or positive-sign signal
        arrvec_t<arr_3d_t> &GC_m,
        const arr_3d_t &psi,
        const arr_3d_t &beta_up,
        const arr_3d_t &beta_dn,
        const arrvec_t<arr_3d_t> &GC_corr,
        const arr_3d_t &G,
        const rng_t &ir,
        const rng_t &jr,
        const rng_t &kr,
        typename std::enable_if<opts::isset(opts, opts::iga) || !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            for (int k = kr.first(); k <= kr.last(); ++k)
            {
              GC_m[d]( pi<d>(i+h, j, k) ) =
              GC_corr[d]( pi<d>(i+h, j, k) ) * where<ix_t>(
                // if
                GC_corr[d]( pi<d>(i+h, j, k) ) > 0,
                // then
                min<ix_t>(1,
                  beta_dn(pi<d>(i,     j, k)),
                  beta_up(pi<d>(i + 1, j, k))
                ),
                // else
                min<ix_t>(1,
                  beta_up(pi<d>(i,     j, k)),
                  beta_dn(pi<d>(i + 1, j, k))
                )
              );
            }
          }
        }
      }
    } // namespace mpdata_fct
  } // namespace formulae
} // namespcae libmpdataxx
