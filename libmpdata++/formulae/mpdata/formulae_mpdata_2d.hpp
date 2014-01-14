/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_common.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_common_2d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_hot_2d.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>

namespace libmpdataxx 
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      // for definition of A and B see eq. 17 in @copybrief Smolarkiewicz_and_Margolin_1998
      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(  // positive sign scalar version
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::pds)>::type* = 0
      ) return_macro(,
	frac<opts>(
	  psi(pi<d>(i+1, j)) - psi(pi<d>(i, j))
	  ,// --------------------------------
	  psi(pi<d>(i+1, j)) + psi(pi<d>(i, j))
	)
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(  // variable-sign scalar version
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::pds)>::type* = 0
      ) return_macro(,
	frac<opts>(
	  abs(psi(pi<d>(i+1, j))) - abs(psi(pi<d>(i, j)))
	  ,// -------------------------------------------
	  abs(psi(pi<d>(i+1, j))) + abs(psi(pi<d>(i, j)))
	) 
      ) 

      template<opts_t opts, int d, class arr_2d_t>
      inline auto A(  // inf. gauge option
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
        (
	  psi(pi<d>(i+1, j)) - psi(pi<d>(i, j))
        ) / ( //---------------------
	  1 + 1
        )
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B( // positive sign signal
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::pds)>::type* = 0
      ) return_macro(,
	frac<opts>( 
	  psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) - psi(pi<d>(i+1, j-1)) - psi(pi<d>(i, j-1))
	  ,// --------------------------------------------------------------------------------
	  psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) + psi(pi<d>(i+1, j-1)) + psi(pi<d>(i, j-1))
	) / 2
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B( // variable-sign signal
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::pds)>::type* = 0
      ) return_macro(,
	frac<opts>( 
	  abs(psi(pi<d>(i+1, j+1))) + abs(psi(pi<d>(i, j+1))) - abs(psi(pi<d>(i+1, j-1))) - abs(psi(pi<d>(i, j-1)))
	  ,// ----------------------------------------------------------------------------------------------------
	  abs(psi(pi<d>(i+1, j+1))) + abs(psi(pi<d>(i, j+1))) + abs(psi(pi<d>(i+1, j-1))) + abs(psi(pi<d>(i, j-1)))
	) / 2
      )

      template<opts_t opts, int d, class arr_2d_t>
      inline auto B( // inf. gauge
        const arr_2d_t &psi, 
        const rng_t &i, 
        const rng_t &j,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0
      ) return_macro(,
        (
	  psi(pi<d>(i+1, j+1)) + psi(pi<d>(i, j+1)) - psi(pi<d>(i+1, j-1)) - psi(pi<d>(i, j-1))
	) / (  // --------------------------------------------------------------------------------
	  1 + 1 + 1 +1
        ) / 2
      )

      //eq. (29a) from @copybrief Smolarkiewicz_and_Margolin_1998
      template <opts_t opts, int dim, class arr_2d_t>
      inline auto antidiff(
        const arr_2d_t &psi, 
        const arrvec_t<arr_2d_t> &GC,
        const arr_2d_t &G, 
        const rng_t &i, 
        const rng_t &j
      ) return_macro(,
        // second order terms
        abs(GC[dim](pi<dim>(i+h, j))) 
        * (1 - abs(GC[dim](pi<dim>(i+h, j))))
        / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j) 
        * A<opts BOOST_PP_COMMA() dim>(psi, i, j) 
        - 
        GC[dim](pi<dim>(i+h, j)) 
        * GC_bar<dim>(GC[dim-1], i, j)
        / G_at_half<opts BOOST_PP_COMMA() dim>(G, i, j) 
        * B<opts BOOST_PP_COMMA() dim>(psi, i, j)
        // third order terms
        + HOT<opts BOOST_PP_COMMA() dim>(psi, GC, G, i, j)
      ) 
    }; // namespace mpdata
  }; // namespace formulae
}; // namespcae libmpdataxx 
