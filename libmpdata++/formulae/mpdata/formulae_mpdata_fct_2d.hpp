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
      template <opts_t opts, class arr_2d_t, class ix_t>
      inline auto beta_up_nominator(
        const arr_2d_t &psi,
        const arr_2d_t &psi_max,
        const arr_2d_t &G,
        const ix_t &i,  
        const ix_t &j
      )
      {
        using std::max;
        return return_helper<ix_t>(
          (  
            max(max(max(max(max(psi_max(i, j), 
                                    psi(i, j+1)),
                      psi(i-1, j)), psi(i, j  )), psi(i+1, j)),
                                    psi(i, j-1)
             ) - psi(i, j)
          ) * formulae::G<opts BOOST_PP_COMMA() 0>(G, i, j) //to make beta up dimensionless when transporting mixing ratios with momentum
        );
      }

      template <opts_t opts, class arr_2d_t, class flx_t>
      inline void beta_up(
        arr_2d_t &b,
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const flx_t &flx,
        const arr_2d_t &G,
        const rng_t &ir,  
        const rng_t &jr
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            b(i, j) = 
            fct_frac<ix_t>(
              beta_up_nominator<opts>(psi, psi_max, G, i, j)
              , // -----------------------------------------------------------
              ( pospart<opts BOOST_PP_COMMA() ix_t>(flx[0](i-h, j))
              - negpart<opts BOOST_PP_COMMA() ix_t>(flx[0](i+h, j)) )  // additional parenthesis so that we first sum
              +                                                        // fluxes in separate dimensions 
              ( pospart<opts BOOST_PP_COMMA() ix_t>(flx[1](i, j-h))    // could be important for accuracy if one of them
              - negpart<opts BOOST_PP_COMMA() ix_t>(flx[1](i, j+h)) )  // is of different magnitude than the other
            );
          }
        }
      } 

      template <opts_t opts, class arr_2d_t, class ix_t>
      inline auto beta_dn_nominator(
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min,
        const arr_2d_t &G, 
        const ix_t &i,
        const ix_t &j 
      )
      {
        using std::min;
        return return_helper<ix_t>(
          (
            psi(i, j)
            - min(min(min(min(min( psi_min(i,j), 
                                       psi(i, j+1)),
                         psi(i-1, j)), psi(i, j  )), psi(i+1, j)),
                                       psi(i, j-1)
            )
          ) * formulae::G<opts BOOST_PP_COMMA() 0>(G, i, j)  //see beta_up_nominator
        );
      } 

      template <opts_t opts, class arr_2d_t, class flx_t>
      inline auto beta_dn(
        arr_2d_t &b, 
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const flx_t &flx,
        const arr_2d_t &G, 
        const rng_t &ir,
        const rng_t &jr
      )
      {
        using ix_t = int;
        for (int i = ir.first(); i <= ir.last(); ++i)
        {
          for (int j = jr.first(); j <= jr.last(); ++j)
          {
            b(i, j) = 
            fct_frac<ix_t>(
              beta_dn_nominator<opts>(psi, psi_min, G, i, j)
              , // ---------------------------------------------------------
              ( pospart<opts BOOST_PP_COMMA() ix_t>(flx[0](i+h, j))
              - negpart<opts BOOST_PP_COMMA() ix_t>(flx[0](i-h, j)) )  //see note in positive sign beta up
              +
              ( pospart<opts BOOST_PP_COMMA() ix_t>(flx[1](i, j+h))
              - negpart<opts BOOST_PP_COMMA() ix_t>(flx[1](i, j-h)) )
            );
          }
        }
      } 

      template <opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto GC_mono( //for variable-sign signal and no infinite gauge option
        const arr_2d_t &psi,
        const arr_2d_t &beta_up,
        const arr_2d_t &beta_dn,
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          GC_corr[d]( pi<d>(i+h, j) ) * where(
            // if
            GC_corr[d]( pi<d>(i+h, j) ) > 0,
            // then
            where(
              // if
              psi(pi<d>(i, j)) > 0,
              // then
              min(1, min(
                beta_dn(pi<d>(i,     j)),
                beta_up(pi<d>(i + 1, j)) 
              )),
              // else
              min(1, min(
                beta_up(pi<d>(i,     j)),
                beta_dn(pi<d>(i + 1, j))
              ))
            ),
            // else
            where(
              // if
              psi(pi<d>(i+1, j)) > 0,
              // then
              min(1, min(
                beta_up(pi<d>(i,     j)),
                beta_dn(pi<d>(i + 1, j))
              )),
              // else
              min(1, min(
                beta_dn(pi<d>(i,     j)),
                beta_up(pi<d>(i + 1, j))
              ))
            )
          )
        );
      } 

      template <opts_t opts, int d, class arr_2d_t, class ix_t>
      inline auto GC_mono( //for infinite gauge option or positive-sign signal
        const arr_2d_t &psi,
        const arr_2d_t &beta_up,
        const arr_2d_t &beta_dn,
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G,
        const ix_t &i,
        const ix_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga) || !opts::isset(opts, opts::abs)>::type* = 0
      )
      {
        return return_helper<ix_t>(
          GC_corr[d]( pi<d>(i+h, j) ) * where(
            // if
            GC_corr[d]( pi<d>(i+h, j) ) > 0, 
            // then
            min(1, min(
              beta_dn(pi<d>(i,     j)),
              beta_up(pi<d>(i + 1, j))
            )),
            // else
            min(1, min(
              beta_up(pi<d>(i,     j)),
              beta_dn(pi<d>(i + 1, j))
            ))
          )
        );
      }
    } // namespace mpdata_fct
  } // namespace formulae
} // namespcae libmpdataxx
