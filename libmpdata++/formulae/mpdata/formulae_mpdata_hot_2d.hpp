/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_common_2d.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {

      //3rd order terms, see  eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998) 
      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_1_helper(
        const arr_2d_t &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
        (
          3 * GC(pi<dim>(i+h, j)) * abs(GC(pi<dim>(i+h, j))) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j)
          - 2 * pow(GC(pi<dim>(i+h, j)), 3) / pow(G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j), 2)
          - GC(pi<dim>(i+h, j))
        ) / 3
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_2_helper(
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j
      ) return_macro(,
        (abs(GC[dim](pi<dim>(i+h, j))) - 2 * pow(GC[dim](pi<dim>(i+h, j)), 2) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j)) 
        / 4 * GC_bar<dim>(GC[dim-1], i, j) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j)
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_1( // positive sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+2, j)) - psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j))
            ,//-----------------------------------------------------------------------------------
            psi(pi<dim>(i+2, j)) + psi(pi<dim>(i+1, j)) + psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_1( // variable sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+2, j))) - abs(psi(pi<dim>(i+1, j))) - abs(psi(pi<dim>(i, j))) + abs(psi(pi<dim>(i-1, j)))
            ,//-------------------------------------------------------------------------------------------------------
            abs(psi(pi<dim>(i+2, j))) + abs(psi(pi<dim>(i+1, j))) + abs(psi(pi<dim>(i, j))) + abs(psi(pi<dim>(i-1, j)))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_1( // inf. gauge option
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs are mutually exclusive");
        ,
	( psi(pi<dim>(i+2, j)) - psi(pi<dim>(i+1, j)) - psi(pi<dim>(i, j)) + psi(pi<dim>(i-1, j)) )
	/ //--------------------------------------------------------------------------------------
	(1 + 1 + 1 + 1)
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_2( // positive sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+1, j+1)) - psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1))            
            ,//-----------------------------------------------------------------------------------------
            psi(pi<dim>(i+1, j+1)) + psi(pi<dim>(i, j+1)) + psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_2( // variable sign signal
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+1, j+1))) - abs(psi(pi<dim>(i, j+1))) - abs(psi(pi<dim>(i+1, j-1))) + abs(psi(pi<dim>(i, j-1)))            
            ,//-------------------------------------------------------------------------------------------------------------
            abs(psi(pi<dim>(i+1, j+1))) + abs(psi(pi<dim>(i, j+1))) + abs(psi(pi<dim>(i+1, j-1))) + abs(psi(pi<dim>(i, j-1)))
        )
      )

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT_2( // inf. gauge option
        const arr_2d_t &psi,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
	(psi(pi<dim>(i+1, j+1)) - psi(pi<dim>(i, j+1)) - psi(pi<dim>(i+1, j-1)) + psi(pi<dim>(i, j-1)))            
	/ //-------------------------------------------------------------------------------------------
	(1 + 1 + 1 + 1)
      )

      // 3rd order term (see eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998 (with G=1))
      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT( // no HOT
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::tot)>::type* = 0
      ) -> decltype(0)
      {
        return 0; //no Higher Order Terms for second accuracy 
      }

      template<opts_t opts, int dim, class arr_2d_t>
      inline auto HOT( // higher order terms
        const arr_2d_t &psi,
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G,
        const rng_t &i,
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::tot)>::type* = 0
      ) return_macro(,
         HOT_1<opts BOOST_PP_COMMA() dim>(psi, i, j) * HOT_1_helper<opts BOOST_PP_COMMA() dim>(GC[dim], G, i, j) 
         + 
         HOT_2<opts BOOST_PP_COMMA() dim>(psi, i, j) * HOT_2_helper<opts BOOST_PP_COMMA() dim>(GC, G, i, j)
      )
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx 
