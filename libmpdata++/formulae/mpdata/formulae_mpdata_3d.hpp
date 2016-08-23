/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_common_3d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_hot_3d.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx
{
  namespace formulae
  {
    namespace mpdata
    {
      // for definition of A and B see eq. 17 in @copybrief Smolarkiewicz_and_Margolin_1998
      template<opts_t opts, int d, class arr_3d_t>
      inline auto A(  // positive sign scalar version
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
	  frac<opts>(
	    psi(pi<d>(i+1, j, k)) - psi(pi<d>(i, j, k))
	  , //-----------------------------------------
	    psi(pi<d>(i+1, j, k)) + psi(pi<d>(i, j, k))
	  )
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto A(  // variable-sign scalar version
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
	  frac<opts>(
	    abs(psi(pi<d>(i+1, j, k))) - abs(psi(pi<d>(i, j, k)))
	  , //---------------------------------------------------
	    abs(psi(pi<d>(i+1, j, k))) + abs(psi(pi<d>(i, j, k)))
	  )
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto A(  // inf. gauge option
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
	(
	  psi(pi<d>(i+1, j, k)) - psi(pi<d>(i, j, k))
	)
	/ //-----------------------------------------
	(
	  1 + 1
	)
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto B1( // positive sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
	  frac<opts>(
	    psi(pi<d>(i+1, j+1, k)) + psi(pi<d>(i, j+1, k)) - psi(pi<d>(i+1, j-1, k)) - psi(pi<d>(i, j-1, k))
	  , //-----------------------------------------------------------------------------------------------
	    psi(pi<d>(i+1, j+1, k)) + psi(pi<d>(i, j+1, k)) + psi(pi<d>(i+1, j-1, k)) + psi(pi<d>(i, j-1, k))
	  ) / 2
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto B1( // variable-sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
	  frac<opts>(
	    abs(psi(pi<d>(i+1, j+1, k)))
          + abs(psi(pi<d>(i, j+1, k)))
          - abs(psi(pi<d>(i+1, j-1, k)))
          - abs(psi(pi<d>(i, j-1, k)))
	  , //--------------------------
	    abs(psi(pi<d>(i+1, j+1, k)))
          + abs(psi(pi<d>(i, j+1, k)))
          + abs(psi(pi<d>(i+1, j-1, k)))
          + abs(psi(pi<d>(i, j-1, k)))
	  ) / 2
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto B1( // inf. gauge
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
	(
	  psi(pi<d>(i+1, j+1, k)) + psi(pi<d>(i, j+1, k)) - psi(pi<d>(i+1, j-1, k)) - psi(pi<d>(i, j-1, k))
	)
	/ //-----------------------------------------------------------------------------------------------
	(
	  1 + 1 + 1 +1
	) / 2
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto B2( // positive sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
	  frac<opts>(
	    psi(pi<d>(i+1, j, k+1)) + psi(pi<d>(i, j, k+1)) - psi(pi<d>(i+1, j, k-1)) - psi(pi<d>(i, j, k-1))
	  , //-----------------------------------------------------------------------------------------------
	    psi(pi<d>(i+1, j, k+1)) + psi(pi<d>(i, j, k+1)) + psi(pi<d>(i+1, j, k-1)) + psi(pi<d>(i, j, k-1))
	  ) / 2
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto B2( // variable-sign signal
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0
      ) return_macro(,
	  frac<opts>(
	    abs(psi(pi<d>(i+1, j, k+1)))
          + abs(psi(pi<d>(i, j, k+1)))
          - abs(psi(pi<d>(i+1, j, k-1)))
          - abs(psi(pi<d>(i, j, k-1)))
	  , //--------------------------
	    abs(psi(pi<d>(i+1, j, k+1)))
          + abs(psi(pi<d>(i, j, k+1)))
          + abs(psi(pi<d>(i+1, j, k-1)))
          + abs(psi(pi<d>(i, j, k-1)))
	  ) / 2
      )

      template<opts_t opts, int d, class arr_3d_t>
      inline auto B2( // inf. gauge
        const arr_3d_t &psi,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "iga & abs options are mutually exclusive");
        ,
	(
	  psi(pi<d>(i+1, j, k+1)) + psi(pi<d>(i, j, k+1)) - psi(pi<d>(i+1, j, k-1)) - psi(pi<d>(i, j, k-1))
	)
	/ //-----------------------------------------------------------------------------------------------
	(
	  1 + 1 + 1 + 1
	) / 2
      )
     
      template <opts_t opts, int dim, class arr_3d_t>
      inline auto antidiff(
        const arr_3d_t &psi,
        const arrvec_t<arr_3d_t> &GC,
        const arr_3d_t &G,
        const rng_t &i,
        const rng_t &j,
        const rng_t &k
      ) return_macro(,
          // second order terms
          abs(GC[dim](pi<dim>(i+h, j, k)))
        * (1 - abs(GC[dim](pi<dim>(i+h, j, k))) / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k))
        * A<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
        - GC[dim](pi<dim>(i+h, j, k))
        * (
            GC_bar1<dim>(GC[dim+1], i, j, k)
          / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k)
          * B1<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
          + GC_bar2<dim>(GC[dim-1], i, j, k)
          / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j, k)
          * B2<opts BOOST_PP_COMMA() dim>(psi, i, j, k)
          )
          // third order terms
        + HOT<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j, k)
      )
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
