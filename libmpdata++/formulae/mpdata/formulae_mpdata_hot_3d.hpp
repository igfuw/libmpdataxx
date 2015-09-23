/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_common_3d.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_1_helper(
        const arr_3d_t &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
          (
            3 * GC(pi<dim>(i+h, j, k)) * abs(GC(pi<dim>(i+h, j, k))) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k)
          - 2 * pow(GC(pi<dim>(i+h, j, k)), 3) / pow(G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k), 2)
          - GC(pi<dim>(i+h, j, k))
          ) / 3
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_2_helper(
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
          (
            abs(GC[dim](pi<dim>(i+h, j, k)))
          - 2 * pow(GC[dim](pi<dim>(i+h, j, k)), 2) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k)
          ) * GC_bar1<dim>(GC[dim+1], i, j, k) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k)
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_3_helper(
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
          (
            abs(GC[dim](pi<dim>(i+h, j, k)))
          - 2 * pow(GC[dim](pi<dim>(i+h, j, k)), 2) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k)
          ) * GC_bar2<dim>(GC[dim-1], i, j, k) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k)
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_4_helper(
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
          -2 * GC[dim](pi<dim>(i+h, j, k)) * GC_bar1<dim>(GC[dim+1], i, j, k) * GC_bar2<dim>(GC[dim-1], i, j, k) / 3
             / pow(G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k), 2)
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_1( // positive sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+2, j, k)) - psi(pi<dim>(i+1, j, k)) - psi(pi<dim>(i, j, k)) + psi(pi<dim>(i-1, j, k))
          , //-------------------------------------------------------------------------------------------------
            psi(pi<dim>(i+2, j, k)) + psi(pi<dim>(i+1, j, k)) + psi(pi<dim>(i, j, k)) + psi(pi<dim>(i-1, j, k))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_1( // variable sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+2, j, k)))
          - abs(psi(pi<dim>(i+1, j, k)))
          - abs(psi(pi<dim>(i, j, k)))
          + abs(psi(pi<dim>(i-1, j, k)))
          , //--------------------------
            abs(psi(pi<dim>(i+2, j, k)))
          + abs(psi(pi<dim>(i+1, j, k)))
          + abs(psi(pi<dim>(i, j, k)))
          + abs(psi(pi<dim>(i-1, j, k)))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_1( // inf. gauge option
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
          static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
          ,
          (
            psi(pi<dim>(i+2, j, k)) - psi(pi<dim>(i+1, j, k)) - psi(pi<dim>(i, j, k)) + psi(pi<dim>(i-1, j, k))
          )
          / //-------------------------------------------------------------------------------------------------
          (
            1 + 1 + 1 + 1
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_2( // positive sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+1, j+1, k)) - psi(pi<dim>(i, j+1, k)) - psi(pi<dim>(i+1, j-1, k)) + psi(pi<dim>(i, j-1, k))
          , //-------------------------------------------------------------------------------------------------------
            psi(pi<dim>(i+1, j+1, k)) + psi(pi<dim>(i, j+1, k)) + psi(pi<dim>(i+1, j-1, k)) + psi(pi<dim>(i, j-1, k))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_2( // variable sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+1, j+1, k)))
          - abs(psi(pi<dim>(i, j+1, k)))
          - abs(psi(pi<dim>(i+1, j-1, k)))
          + abs(psi(pi<dim>(i, j-1, k)))           
          , //----------------------------
            abs(psi(pi<dim>(i+1, j+1, k)))
          + abs(psi(pi<dim>(i, j+1, k)))
          + abs(psi(pi<dim>(i+1, j-1, k)))
          + abs(psi(pi<dim>(i, j-1, k)))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_2( // inf. gauge option
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
          static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
          ,
          (
            psi(pi<dim>(i+1, j+1, k))
          - psi(pi<dim>(i, j+1, k))
          - psi(pi<dim>(i+1, j-1, k))
          + psi(pi<dim>(i, j-1, k))
          )
          / //-----------------------
          (
            1 + 1 + 1 + 1
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_3( // positive sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+1, j, k+1))
          - psi(pi<dim>(i, j, k+1))
          - psi(pi<dim>(i+1, j, k-1))
          + psi(pi<dim>(i, j, k-1))           
          , //-----------------------
            psi(pi<dim>(i+1, j, k+1))
          + psi(pi<dim>(i, j, k+1))
          + psi(pi<dim>(i+1, j, k-1))
          + psi(pi<dim>(i, j, k-1))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_3( // variable sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+1, j, k+1)))
          - abs(psi(pi<dim>(i, j, k+1)))
          - abs(psi(pi<dim>(i+1, j, k-1)))
          + abs(psi(pi<dim>(i, j, k-1)))           
          , //----------------------------
            abs(psi(pi<dim>(i+1, j, k+1)))
          + abs(psi(pi<dim>(i, j, k+1)))
          + abs(psi(pi<dim>(i+1, j, k-1)))
          + abs(psi(pi<dim>(i, j, k-1)))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_3( // inf. gauge option
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
          static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
          ,
          (
            psi(pi<dim>(i+1, j, k+1))
          - psi(pi<dim>(i, j, k+1))
          - psi(pi<dim>(i+1, j, k-1))
          + psi(pi<dim>(i, j, k-1))
          )           
          / //-----------------------
          (
            1 + 1 + 1 + 1
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_4( // positive sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            psi(pi<dim>(i+1, j+1, k+1))
          + psi(pi<dim>(i+1, j-1, k-1))
          - psi(pi<dim>(i+1, j+1, k-1))
          - psi(pi<dim>(i+1, j-1, k+1))           
          + psi(pi<dim>(i, j+1, k+1))
          + psi(pi<dim>(i, j-1, k-1))
          - psi(pi<dim>(i, j+1, k-1))           
          - psi(pi<dim>(i, j-1, k+1))
          , //-------------------------
            psi(pi<dim>(i+1, j+1, k+1))
          + psi(pi<dim>(i+1, j-1, k-1))
          + psi(pi<dim>(i+1, j+1, k-1))
          + psi(pi<dim>(i+1, j-1, k+1))           
          + psi(pi<dim>(i, j+1, k+1))
          + psi(pi<dim>(i, j-1, k-1))
          + psi(pi<dim>(i, j+1, k-1))           
          + psi(pi<dim>(i, j-1, k+1))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_4( // variable sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
          frac<opts>(
            abs(psi(pi<dim>(i+1, j+1, k+1)))
          + abs(psi(pi<dim>(i+1, j-1, k-1)))
          - abs(psi(pi<dim>(i+1, j+1, k-1)))
          - abs(psi(pi<dim>(i+1, j-1, k+1)))           
          + abs(psi(pi<dim>(i, j+1, k+1)))
          + abs(psi(pi<dim>(i, j-1, k-1)))
          - abs(psi(pi<dim>(i, j+1, k-1)))          
          - abs(psi(pi<dim>(i, j-1, k+1)))
          , //------------------------------
            abs(psi(pi<dim>(i+1, j+1, k+1)))
          + abs(psi(pi<dim>(i+1, j-1, k-1)))
          + abs(psi(pi<dim>(i+1, j+1, k-1)))
          + abs(psi(pi<dim>(i+1, j-1, k+1)))           
          + abs(psi(pi<dim>(i, j+1, k+1)))
          + abs(psi(pi<dim>(i, j-1, k-1)))
          + abs(psi(pi<dim>(i, j+1, k-1)))          
          + abs(psi(pi<dim>(i, j-1, k+1)))
          )
      )

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT_4( // inf. gauge option
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
          static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
          ,
          (
            psi(pi<dim>(i+1, j+1, k+1))
          + psi(pi<dim>(i+1, j-1, k-1))
          - psi(pi<dim>(i+1, j+1, k-1))
          - psi(pi<dim>(i+1, j-1, k+1))
          + psi(pi<dim>(i, j+1, k+1))
          + psi(pi<dim>(i, j-1, k-1))           
          - psi(pi<dim>(i, j+1, k-1))
          - psi(pi<dim>(i, j-1, k+1))
          )
          / //---------------------------
          (
            1 + 1 + 1 + 1 + 1 + 1 + 1 + 1
          )
      )

      // 3rd order term (see eq. (36) from @copybrief Smolarkiewicz_and_Margolin_1998 (with G=1))
      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT( // no HOT
        const arr_3d_t &psi,
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::tot)>::type* = 0
      ) -> decltype(0)
      {
        return 0; //no Higher Order Terms for second accuracy
      }

      template<opts_t opts, int dim, class arr_3d_t>
      inline auto HOT( // higher order terms
        const arr_3d_t &psi,
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::tot)>::type* = 0
      ) return_macro(,
          HOT_1<opts BOOST_PP_COMMA() dim>(psi, i, j, k) * HOT_1_helper<opts BOOST_PP_COMMA() dim>(GC[dim], G, i, j, k)
        + HOT_2<opts BOOST_PP_COMMA() dim>(psi, i, j, k) * HOT_2_helper<opts BOOST_PP_COMMA() dim>(GC, G, i, j, k)
        + HOT_3<opts BOOST_PP_COMMA() dim>(psi, i, j, k) * HOT_3_helper<opts BOOST_PP_COMMA() dim>(GC, G, i, j, k)
        + HOT_4<opts BOOST_PP_COMMA() dim>(psi, i, j, k) * HOT_4_helper<opts BOOST_PP_COMMA() dim>(GC, G, i, j, k)
      )


    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
