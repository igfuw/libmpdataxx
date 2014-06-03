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
      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_up_nominator(
        const arr_2d_t &psi,
        const arr_2d_t &psi_max,
        const arr_2d_t &G,
        const rng_t i,  
        const rng_t j
      ) return_macro(,
          (  
            max(max(max(max(max(psi_max(pi<dim>(i, j)), 
                                    psi(pi<dim>(i, j+1))),
             psi(pi<dim>(i-1, j))), psi(pi<dim>(i, j  ))), psi(pi<dim>(i+1, j))),
                                    psi(pi<dim>(i, j-1))
             ) - psi(pi<dim>(i, j))
          ) * formulae::G<opts BOOST_PP_COMMA() dim>(G, i, j) //to make beta up dimensionless when transporting mixing ratios with momentum
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_up( //positive sign signal
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G,
        const rng_t i,  
        const rng_t j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        fct_frac(
          beta_up_nominator<opts, dim>(psi, psi_max, G, i, j)
          , // -----------------------------------------------------------
          ( pospart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) * psi(pi<dim>(i-1, j)) 
          - negpart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) * psi(pi<dim>(i+1, j)) )  // additional parenthesis so that we first sum
          +                                                                          // fluxes in separate dimensions 
          ( pospart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) * psi(pi<dim>(i, j-1))    // could be important for accuracy if one of them
          - negpart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) * psi(pi<dim>(i, j+1)) )  // is of different magnitude than the other
        ) 
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_up( //variable-sign signal
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G,
        const rng_t i,  
        const rng_t j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        fct_frac(
          beta_up_nominator<opts, dim>(psi, psi_max, G, i, j)
          , // --------------------------------------------------------------------------
          ( pospart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) * pospart<opts>(psi(pi<dim>(i-1, j))) 
          - negpart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) * pospart<opts>(psi(pi<dim>(i+1, j)))
          - pospart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) * negpart<opts>(psi(pi<dim>(i,   j)))
          + negpart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) * negpart<opts>(psi(pi<dim>(i,   j))) ) // see note in positive sign beta up
          +
          ( pospart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) * pospart<opts>(psi(pi<dim>(i, j-1)))
          - negpart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) * pospart<opts>(psi(pi<dim>(i, j+1)))
          - pospart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) * negpart<opts>(psi(pi<dim>(i, j  )))
          + negpart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) * negpart<opts>(psi(pi<dim>(i, j  ))) )
        ) 
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_up( //inf. gauge option
        const arr_2d_t &psi,
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G,
        const rng_t i,  
        const rng_t j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
        fct_frac(
          beta_up_nominator<opts, dim>(psi, psi_max, G, i, j)
          , // -------------------------------------------
          ( pospart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) /* *1 */ 
          - negpart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) /* *1 */) // see note in positive sign beta up
          +
          ( pospart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) /* *1 */
          - negpart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) /* *1 */)
        ) 
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_dn_nominator(
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min,
        const arr_2d_t &G, 
        const rng_t i,
        const rng_t j 
      ) return_macro(,
          (
            psi(pi<dim>(i, j))
            - min(min(min(min(min( psi_min(pi<dim>(i,j)), 
                                       psi(pi<dim>(i, j+1))),
                psi(pi<dim>(i-1, j))), psi(pi<dim>(i, j  ))), psi(pi<dim>(i+1, j))),
                                       psi(pi<dim>(i, j-1))
            )
          ) * formulae::G<opts BOOST_PP_COMMA() dim>(G, i, j)  //see beta_up_nominator
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_dn( //positive-sign signal
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G, 
        const rng_t i,
        const rng_t j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
        fct_frac(
          beta_dn_nominator<opts, dim>(psi, psi_min, G, i, j)
	  , // ---------------------------------------------------------
          ( pospart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) * psi(pi<dim>(i, j))
          - negpart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) * psi(pi<dim>(i, j)) )  //see note in positive sign beta up
          +
          ( pospart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) * psi(pi<dim>(i, j))
          - negpart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) * psi(pi<dim>(i, j)) )
        ) 
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_dn( //variable-sign signal
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G, 
        const rng_t i,
        const rng_t j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
        fct_frac(
          beta_dn_nominator<opts, dim>(psi, psi_min, G, i, j)
          , // --------------------------
          ( pospart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) * pospart<opts>(psi(pi<dim>(i,   j)))
          - negpart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) * pospart<opts>(psi(pi<dim>(i,   j)))
          - pospart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) * negpart<opts>(psi(pi<dim>(i-1, j)))
          + negpart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) * negpart<opts>(psi(pi<dim>(i+1, j))) )  //see note in positive sign beta up
          +
          ( pospart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) * pospart<opts>(psi(pi<dim>(i,   j)))
          - negpart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) * pospart<opts>(psi(pi<dim>(i,   j)))
          - pospart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) * negpart<opts>(psi(pi<dim>(i, j-1)))
          + negpart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) * negpart<opts>(psi(pi<dim>(i, j+1))) )
        ) 
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto beta_dn( //inf. gauge option
        const arr_2d_t &psi, 
        const arr_2d_t &psi_min, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G, 
        const rng_t i,
        const rng_t j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0 
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
        fct_frac(
          beta_dn_nominator<opts, dim>(psi, psi_min, G, i, j)
          , // -------------------------------------------
          ( pospart<opts>(GC_corr[dim-0](pi<dim>(i+h, j))) /* *1 */
          - negpart<opts>(GC_corr[dim-0](pi<dim>(i-h, j))) /* *1 */)  //see note in positive sign beta up
          +
          ( pospart<opts>(GC_corr[dim-1](pi<dim>(i, j+h))) /* *1 */
          - negpart<opts>(GC_corr[dim-1](pi<dim>(i, j-h))) /* *1 */)
        ) 
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto GC_mono( //for variable-sign signal and no infinite gauge option
        const arr_2d_t &psi,
        const arr_2d_t &psi_min, // from before the first iteration
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G,
        const rng_t i,
        const rng_t j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        GC_corr[dim]( pi<dim>(i+h, j) ) * where( // TODO: is it possible to implement it without where()?
          // if
          GC_corr[dim]( pi<dim>(i+h, j) ) > 0, // >= ?
          // then
          where(
            // if
            psi(pi<dim>(i, j)) > 0, // TODO: wouldn't it be better with >= ?
            // then
	    min(1, min(
              beta_dn<opts, dim>(psi, psi_min, GC_corr, G, i,     j),
              beta_up<opts, dim>(psi, psi_max, GC_corr, G, i + 1, j) 
	    )),
            // else
	    min(1, min(
	      beta_up<opts, dim>(psi, psi_max, GC_corr, G, i,     j),
	      beta_dn<opts, dim>(psi, psi_min, GC_corr, G, i + 1, j)
	    ))
          ),
          // else
          where(
            // if
            psi(pi<dim>(i+1, j)) > 0, // TODO: what if crossing zero?
            // then
	    min(1, min(
	      beta_up<opts, dim>(psi, psi_max, GC_corr, G, i,     j),
	      beta_dn<opts, dim>(psi, psi_min, GC_corr, G, i + 1, j)
	    )),
            // else
	    min(1, min(
	      beta_dn<opts, dim>(psi, psi_min, GC_corr, G, i,     j),
	      beta_up<opts, dim>(psi, psi_max, GC_corr, G, i + 1, j)
	    ))
          )
        )
      ) 

      template <opts_t opts, int dim, class arr_2d_t>
      inline auto GC_mono( //for infinite gauge option or positive-sign signal
        const arr_2d_t &psi,
        const arr_2d_t &psi_min, // from before the first iteration
        const arr_2d_t &psi_max, // from before the first iteration
        const arrvec_t<arr_2d_t> &GC_corr,
        const arr_2d_t &G,
        const rng_t i,
        const rng_t j,
        typename std::enable_if<opts::isset(opts, opts::iga) || !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
        GC_corr[dim]( pi<dim>(i+h, j) ) * where(
          // if
          GC_corr[dim]( pi<dim>(i+h, j) ) >= 0, // TODO: what about where(C!=0, where()...) - could be faster for fields with a lot of zeros? (as an option?)
          // then
	  min(1, min(
	    beta_dn<opts, dim>(psi, psi_min, GC_corr, G, i,     j),
	    beta_up<opts, dim>(psi, psi_max, GC_corr, G, i + 1, j)
	  )),
          // else
	  min(1, min(
	    beta_up<opts, dim>(psi, psi_max, GC_corr, G, i,     j),
	    beta_dn<opts, dim>(psi, psi_min, GC_corr, G, i + 1, j)
	  ))
        )
      )
    }; // namespace mpdata_fct
  }; // namespace formulae
}; // namespcae libmpdataxx
