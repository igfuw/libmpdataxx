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

      template <int d, class arr_2d_t>
      inline auto beta_up_nominator(
        const arr_2d_t &psi,
        const arr_2d_t &psi_max,
        const rng_t i,  
        const rng_t j
      ) return_macro(,
          max(max(max(max(max(psi_max(pi<d>(i, j)), 
                                 psi(pi<d>(i, j+1))),
            psi(pi<d>(i-1, j))), psi(pi<d>(i, j  ))), psi(pi<d>(i+1, j))),
                                 psi(pi<d>(i, j-1))
          ) 
          - psi(pi<d>(i, j))
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_up( //positive sign signal
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,  
        const rng_t j,
        typename std::enable_if<opt_set(opts, pds)>::type* = 0
      ) return_macro(,
        frac<opts>(
          beta_up_nominator<d>(psi, psi_max, i, j)
          ,// -----------------------------------------------------------
            pospart<opts>(C_corr[d-0](pi<d>(i-h, j))) * psi(pi<d>(i-1, j)) 
          - negpart<opts>(C_corr[d-0](pi<d>(i+h, j))) * psi(pi<d>(i+1, j))

          + pospart<opts>(C_corr[d-1](pi<d>(i, j-h))) * psi(pi<d>(i, j-1))
          - negpart<opts>(C_corr[d-1](pi<d>(i, j+h))) * psi(pi<d>(i, j+1))
        ) 
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_up( //variable-sign signal
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,  
        const rng_t j,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0
      ) return_macro(,
        frac<opts>(
          beta_up_nominator<d>(psi, psi_max, i, j)
          ,// --------------------------------------------------------------------------
            pospart<opts>(C_corr[d-0](pi<d>(i-h, j))) * pospart<opts>(psi(pi<d>(i-1, j))) 
          - negpart<opts>(C_corr[d-0](pi<d>(i+h, j))) * pospart<opts>(psi(pi<d>(i+1, j)))
          - pospart<opts>(C_corr[d-0](pi<d>(i+h, j))) * negpart<opts>(psi(pi<d>(i,   j)))
          + negpart<opts>(C_corr[d-0](pi<d>(i-h, j))) * negpart<opts>(psi(pi<d>(i,   j)))

          + pospart<opts>(C_corr[d-1](pi<d>(i, j-h))) * pospart<opts>(psi(pi<d>(i, j-1)))
          - negpart<opts>(C_corr[d-1](pi<d>(i, j+h))) * pospart<opts>(psi(pi<d>(i, j+1)))
          - pospart<opts>(C_corr[d-1](pi<d>(i, j+h))) * negpart<opts>(psi(pi<d>(i, j  )))
          + negpart<opts>(C_corr[d-1](pi<d>(i, j-h))) * negpart<opts>(psi(pi<d>(i, j  )))
        ) 
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_up( //inf. gauge option
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,  
        const rng_t j,
        typename std::enable_if<opt_set(opts, iga)>::type* = 0
      ) return_macro(,
        frac<opts>(
          beta_up_nominator<d>(psi, psi_max, i, j)
          ,// -------------------------------------------
            pospart<opts>(C_corr[d-0](pi<d>(i-h, j))) * 1 
          - negpart<opts>(C_corr[d-0](pi<d>(i+h, j))) * 1

          + pospart<opts>(C_corr[d-1](pi<d>(i, j-h))) * 1
          - negpart<opts>(C_corr[d-1](pi<d>(i, j+h))) * 1
        ) 
      ) 

      template <int d, class arr_2d_t>
      inline auto beta_dn_nominator(
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min,
        const rng_t i,
        const rng_t j 
      ) return_macro(,
          psi(pi<d>(i, j))
          - min(min(min(min(min(psi_min(pi<d>(i,j)), 
                                   psi(pi<d>(i, j+1))),
              psi(pi<d>(i-1, j))), psi(pi<d>(i, j  ))), psi(pi<d>(i+1, j))),
                                   psi(pi<d>(i, j-1))
        )
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_dn( //positive-sign signal
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        typename std::enable_if<opt_set(opts, pds)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_dn_nominator<d>(psi, psi_min, i, j)
          ,// ---------------------------------------------------------
            pospart<opts>(C_corr[d-0](pi<d>(i+h, j))) * psi(pi<d>(i, j))
          - negpart<opts>(C_corr[d-0](pi<d>(i-h, j))) * psi(pi<d>(i, j))

          + pospart<opts>(C_corr[d-1](pi<d>(i, j+h))) * psi(pi<d>(i, j))
          - negpart<opts>(C_corr[d-1](pi<d>(i, j-h))) * psi(pi<d>(i, j))
        ) 
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto beta_dn( //variable-sign signal
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_dn_nominator<d>(psi, psi_min, i, j)
          ,// --------------------------
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
      inline auto beta_dn( //inf. gauge option
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        typename std::enable_if<opt_set(opts, iga)>::type* = 0 
      ) return_macro(,
        frac<opts>(
          beta_dn_nominator<d>(psi, psi_min, i, j)
          ,// -------------------------------------------
            pospart<opts>(C_corr[d-0](pi<d>(i+h, j))) * 1
          - negpart<opts>(C_corr[d-0](pi<d>(i-h, j))) * 1

          + pospart<opts>(C_corr[d-1](pi<d>(i, j+h))) * 1
          - negpart<opts>(C_corr[d-1](pi<d>(i, j-h))) * 1
        ) 
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto C_mono( //for variable-sign signal and no infinite gauge option
        const arr_2d_t &psi,
        const arr_2d_t &psi_min, // from before the first iteration
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        typename std::enable_if<!opt_set(opts, iga) && !opt_set(opts, pds)>::type* = 0
      ) return_macro(,
        C_corr[d]( pi<d>(i+h, j) ) * where( // TODO: is it possible to implement it without where()?
          // if
          C_corr[d]( pi<d>(i+h, j) ) > 0, // >= ?
          // then
          where(
            // if
            psi(pi<d>(i, j)) > 0, // TODO: wouldn't it be better with >= ?
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
          ),
          // else
          where(
            // if
            psi(pi<d>(i+1, j)) > 0, // TODO: what if crossing zero?
            // then
	    min(1, min(
	      beta_up<opts, d>(psi, psi_max, C_corr, i,     j),
	      beta_dn<opts, d>(psi, psi_min, C_corr, i + 1, j)
	    )),
            // else
	    min(1, min(
	      beta_dn<opts, d>(psi, psi_min, C_corr, i,     j),
	      beta_up<opts, d>(psi, psi_max, C_corr, i + 1, j)
	    ))
          )
        )
      ) 

      template <opts_t opts, int d, class arr_2d_t>
      inline auto C_mono( //for infinite gauge option or positive-sign signal
        const arr_2d_t &psi,
        const arr_2d_t &psi_min, // from before the first iteration
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &C_corr,
        const rng_t i,
        const rng_t j,
        typename std::enable_if<opt_set(opts, iga) || opt_set(opts, pds)>::type* = 0
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
