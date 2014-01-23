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

      template <int d, class arr_3d_t>
      inline auto beta_up_nominator(
        const arr_3d_t &psi,
        const arr_3d_t &psi_max,
        const rng_t i,  
        const rng_t j,
        const rng_t k
      ) return_macro(,
          max(max(max(max(max(max(max(
                      psi_max(pi<d>(i, j, k)),
                      psi(pi<d>(i, j, k))),
                      psi(pi<d>(i+1, j, k))),
                      psi(pi<d>(i-1, j, k))),
                      psi(pi<d>(i, j+1, k))),
                      psi(pi<d>(i, j-1, k))),
                      psi(pi<d>(i, j, k+1))),
                      psi(pi<d>(i, j, k-1))
          ) 
        - psi(pi<d>(i, j, k))
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto beta_up( //positive sign signal
        const arr_3d_t &psi,
        const arr_3d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,  
        const rng_t j,
        const rng_t k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        frac<opts>(
          beta_up_nominator<d>(psi, psi_max, i, j, k)
        , //----------------------------------------------------------------------
          ( pospart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * psi(pi<d>(i-1, j, k)) 
          - negpart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * psi(pi<d>(i+1, j, k)) )  // additional parenthesis so that we first sum
          +                                                                   // fluxes in separate dimensions 
          ( pospart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * psi(pi<d>(i, j-1, k))    // could be important for accuracy if one of them
          - negpart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * psi(pi<d>(i, j+1, k)) )  // is of different magnitude than the other
          +                                                                   // fluxes in separate dimensions 
          ( pospart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * psi(pi<d>(i, j, k-1))    
          - negpart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * psi(pi<d>(i, j, k+1)) )
        ) 
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto beta_up( //variable-sign signal
        const arr_3d_t &psi,
        const arr_3d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,  
        const rng_t j,
        const rng_t k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        frac<opts>(
          beta_up_nominator<d>(psi, psi_max, i, j, k)
        , //-------------------------------------------------------------------------------------
          ( pospart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * pospart<opts>(psi(pi<d>(i-1, j, k))) 
          - negpart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * pospart<opts>(psi(pi<d>(i+1, j, k)))
          - pospart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * negpart<opts>(psi(pi<d>(i,   j, k)))
          + negpart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * negpart<opts>(psi(pi<d>(i,   j, k))) ) // see note in positive sign beta up
          +
          ( pospart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * pospart<opts>(psi(pi<d>(i, j-1, k)))
          - negpart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * pospart<opts>(psi(pi<d>(i, j+1, k)))
          - pospart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * negpart<opts>(psi(pi<d>(i, j  , k)))
          + negpart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * negpart<opts>(psi(pi<d>(i, j  , k))) )
          +
          ( pospart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * pospart<opts>(psi(pi<d>(i, j, k-1)))
          - negpart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * pospart<opts>(psi(pi<d>(i, j, k+1)))
          - pospart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * negpart<opts>(psi(pi<d>(i, j  , k)))
          + negpart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * negpart<opts>(psi(pi<d>(i, j  , k))) )
        ) 
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto beta_up( //inf. gauge option
        const arr_3d_t &psi,
        const arr_3d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,  
        const rng_t j,
        const rng_t k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
        frac<opts>(
          beta_up_nominator<d>(psi, psi_max, i, j, k)
        , //--------------------------------------------------
          ( pospart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * 1 
          - negpart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * 1 ) // see note in positive sign beta up
          +
          ( pospart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * 1
          - negpart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * 1 )
          +
          ( pospart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * 1
          - negpart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * 1 )
        ) 
      ) 

      template <int d, class arr_3d_t>
      inline auto beta_dn_nominator(
        const arr_3d_t &psi, 
        const arr_3d_t &psi_min,
        const rng_t i,
        const rng_t j,
        const rng_t k
      ) return_macro(,
          psi(pi<d>(i, j, k))
          - min(min(min(min(min(min(min(
                    psi_min(pi<d>(i, j, k)), 
                    psi(pi<d>(i, j+1, k))),
                    psi(pi<d>(i-1, j, k))),
                    psi(pi<d>(i, j, k  ))),
                    psi(pi<d>(i+1, j, k))),
                    psi(pi<d>(i, j, k+1))),
                    psi(pi<d>(i, j, k-1))),
                    psi(pi<d>(i, j-1, k))
        )
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto beta_dn( //positive-sign signal
        const arr_3d_t &psi, 
        const arr_3d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        const rng_t k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_dn_nominator<d>(psi, psi_min, i, j, k)
        , //--------------------------------------------------------------------
          ( pospart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * psi(pi<d>(i, j, k))
          - negpart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * psi(pi<d>(i, j, k)) )  //see note in positive sign beta up
          +
          ( pospart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * psi(pi<d>(i, j, k))
          - negpart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * psi(pi<d>(i, j, k)) )
          +
          ( pospart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * psi(pi<d>(i, j, k))
          - negpart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * psi(pi<d>(i, j, k)) )
        ) 
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto beta_dn( //variable-sign signal
        const arr_3d_t &psi, 
        const arr_3d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        const rng_t k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_dn_nominator<d>(psi, psi_min, i, j, k)
        , //-----------------------------------------------------------------------------------
          ( pospart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * pospart<opts>(psi(pi<d>(i,   j, k)))
          - negpart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * pospart<opts>(psi(pi<d>(i,   j, k)))
          - pospart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * negpart<opts>(psi(pi<d>(i-1, j, k)))
          + negpart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * negpart<opts>(psi(pi<d>(i+1, j, k))) )  //see note in positive sign beta up
          +
          ( pospart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * pospart<opts>(psi(pi<d>(i,   j, k)))
          - negpart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * pospart<opts>(psi(pi<d>(i,   j, k)))
          - pospart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * negpart<opts>(psi(pi<d>(i, j-1, k)))
          + negpart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * negpart<opts>(psi(pi<d>(i, j+1, k))) )
          +
          ( pospart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * pospart<opts>(psi(pi<d>(i,   j, k)))
          - negpart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * pospart<opts>(psi(pi<d>(i,   j, k)))
          - pospart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * negpart<opts>(psi(pi<d>(i, j, k-1)))
          + negpart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * negpart<opts>(psi(pi<d>(i, j, k+1))) )
        ) 
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto beta_dn( //inf. gauge option
        const arr_3d_t &psi, 
        const arr_3d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        const rng_t k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs are mutually exclusive");
        ,
        frac<opts>(
          beta_dn_nominator<d>(psi, psi_min, i, j, k)
        , //--------------------------------------------------
          ( pospart<opts>(C_corr[d+0](pi<d>(i+h, j, k))) * 1
          - negpart<opts>(C_corr[d+0](pi<d>(i-h, j, k))) * 1 )  //see note in positive sign beta up
          +
          ( pospart<opts>(C_corr[d+1](pi<d>(i, j+h, k))) * 1
          - negpart<opts>(C_corr[d+1](pi<d>(i, j-h, k))) * 1 )
          +
          ( pospart<opts>(C_corr[d+2](pi<d>(i, j, k+h))) * 1
          - negpart<opts>(C_corr[d+2](pi<d>(i, j, k-h))) * 1 )
        ) 
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto C_mono( //for variable-sign signal and no infinite gauge option
        const arr_3d_t &psi,
        const arr_3d_t &psi_min, // from before the first iteration
        const arr_3d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        const rng_t k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        C_corr[d]( pi<d>(i+h, j, k) ) * where( // TODO: is it possible to implement it without where()?
          // if
          C_corr[d]( pi<d>(i+h, j, k) ) > 0, // >= ?
          // then
          where(
            // if
            psi(pi<d>(i, j, k)) > 0, // TODO: wouldn't it be better with >= ?
            // then
	    min(1, min(
              beta_dn<opts, d>(psi, psi_min, C_corr, i,     j, k),
              beta_up<opts, d>(psi, psi_max, C_corr, i + 1, j, k) 
	    )),
            // else
	    min(1, min(
	      beta_up<opts, d>(psi, psi_max, C_corr, i,     j, k),
	      beta_dn<opts, d>(psi, psi_min, C_corr, i + 1, j, k)
	    ))
          ),
          // else
          where(
            // if
            psi(pi<d>(i+1, j, k)) > 0, // TODO: what if crossing zero?
            // then
	    min(1, min(
	      beta_up<opts, d>(psi, psi_max, C_corr, i,     j, k),
	      beta_dn<opts, d>(psi, psi_min, C_corr, i + 1, j, k)
	    )),
            // else
	    min(1, min(
	      beta_dn<opts, d>(psi, psi_min, C_corr, i,     j, k),
	      beta_up<opts, d>(psi, psi_max, C_corr, i + 1, j, k)
	    ))
          )
        )
      ) 

      template <opts_t opts, int d, class arr_3d_t>
      inline auto C_mono( //for infinite gauge option or positive-sign signal
        const arr_3d_t &psi,
        const arr_3d_t &psi_min, // from before the first iteration
        const arr_3d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_3d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        const rng_t k,
        typename std::enable_if<opts::isset(opts, opts::iga) || !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        C_corr[d]( pi<d>(i+h, j, k) ) * where(
          // if
          C_corr[d]( pi<d>(i+h, j, k) ) >= 0, // TODO: what about where(C!=0, where()...) - could be faster for fields with a lot of zeros? (as an option?)
          // then
	  min(1, min(
	    beta_dn<opts, d>(psi, psi_min, C_corr, i,     j, k),
	    beta_up<opts, d>(psi, psi_max, C_corr, i + 1, j, k)
	  )),
          // else
	  min(1, min(
	    beta_up<opts, d>(psi, psi_max, C_corr, i,     j, k),
	    beta_dn<opts, d>(psi, psi_min, C_corr, i + 1, j, k)
	  ))
        )
      )
    }; // namespace mpdata_fct
  }; // namespace formulae
}; // namespcae libmpdataxx
