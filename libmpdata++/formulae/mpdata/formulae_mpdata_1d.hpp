/** @file
* @copyright University of Warsaw
* @section LICENSE
* GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
*/

#pragma once

#include <libmpdata++/formulae/mpdata/formulae_mpdata_hot_1d.hpp>
#include <libmpdata++/formulae/mpdata/formulae_mpdata_dfl_1d.hpp>

namespace libmpdataxx
{ 
  namespace formulae 
  { 
    namespace mpdata 
    {
      //eq. (17) from @copybrief Smolarkiewicz_and_Margolin_1998
      template <opts_t opts, class arr_1d_t>
      inline auto A(  // positive-sign signal version (no need for abs())
        const arr_1d_t &psi, 
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && !opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
	frac<opts>( 
   	  psi(i+1) - psi(i)
	  ,// -------------
          psi(i+1) + psi(i)
	)
      )

      template <opts_t opts, class arr_1d_t>
      inline auto A(   // variable-sign signal version (hence the need for abs())
        const arr_1d_t &psi, 
        const rng_t &i,
        typename std::enable_if<!opts::isset(opts, opts::iga) && opts::isset(opts, opts::abs)>::type* = 0 
      ) return_macro(,
	frac<opts>( 
          abs(psi(i+1)) - abs(psi(i))
	  ,// -----------------------
	  abs(psi(i+1)) + abs(psi(i))
	) 
      )

      template <opts_t opts, class arr_1d_t>
      inline auto A( // infinite gauge version (for both signed and variable-sign fields), (sum of psi -> sum of 1)
        const arr_1d_t &psi, 
        const rng_t &i,
        typename std::enable_if<opts::isset(opts, opts::iga)>::type* = 0 // enabled if iga == true
      ) return_macro(
        static_assert(!opts::isset(opts, opts::abs), "abs & iga options are mutually exclusive");
        ,
	(psi(i+1) - psi(i)) 
        / //---------------
        (1 + 1)
      )
      // eq 29a from @copybrief Smolarkiewicz_and_Margolin_1998
      template<opts_t opts, class arr_1d_t>
      inline auto antidiff( // antidiffusive velocity
        const arr_1d_t &psi, 
        const arr_1d_t &GC,
        const arr_1d_t &G,
        const rng_t &i
      ) return_macro(,
        // second-order terms
        abs(GC(i+h)) 
        * (1 - abs(GC(i+h)) / ((formulae::G<opts>(G, i) + formulae::G<opts>(G, i+1)) / 2)) 
        * A<opts>(psi, i) 
        // third-order terms
        + HOT<opts>(psi, GC, G, i) //higher order term
        // divergent flow terms
        + DFL<opts>(psi, GC, G, i) //divergent flow correction
      )
    } // namespace mpdata
  } // namespace formulae
} // namespcae libmpdataxx
